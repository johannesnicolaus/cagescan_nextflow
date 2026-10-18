"""READ1 5' end and spliced-leader (SL) helpers shared by bam_to_ctss.py and assign_pairs_to_tss.py.

Leaders are given as a comma-separated list, optionally named: "SEQ" or "SL1=SEQ1,SL2a=SEQ2,...". Unnamed leaders
are called SL (one leader) or SL1, SL2, ... (several). Pairs whose leader was trimmed before alignment
(SPLIT_SPLICED_LEADER) are renamed "SL.<leader>:<name>" or "SLG.<leader>:<name>" (extra G before the leader);
"SL:" and "SLG:" (one leader, earlier pipeline versions) are recognised too.
"""

import re

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SL_TRIMMED = ("SL.", "SLG.", "SL:", "SLG:")
MAX_SL_CLIP = 3  # longest 5' soft clip of an SL-trimmed READ1 that is added back to its 5' end
MD_FIRST = re.compile(r"^0[ACGTN]")
MD_LAST = re.compile(r"[ACGTN]0$")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_leaders(spec):
    """'SEQ' or 'NAME=SEQ,...' -> list of (name, SEQ); [] if empty."""
    entries = [e.strip() for e in (spec or "").split(",") if e.strip()]
    leaders = []
    for i, entry in enumerate(entries, 1):
        name, _, seq = entry.rpartition("=")
        leaders.append((name or ("SL" if len(entries) == 1 else f"SL{i}"), seq.upper()))
    names = [name for name, _ in leaders]
    if len(set(names)) < len(names):
        raise ValueError(f"duplicate spliced-leader names in {spec!r}: {names}")
    return leaders


def five_prime(read):
    """Return (pos0, strand, softclipped 5' sequence, first aligned bases), all in read orientation."""
    cigar = [op for op in (read.cigartuples or []) if op[0] != 5]  # hard clips are not in the sequence
    seq = read.query_sequence or ""
    if read.is_reverse:
        clip = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
        oriented = revcomp(seq)
        return read.reference_end - 1, "-", oriented[:clip], oriented[clip:clip + 8]
    clip = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    return read.reference_start, "+", seq[:clip], seq[clip:clip + 8]


def shift_downstream(pos, strand, k):
    return pos + k if strand == "+" else pos - k


def trimmed_leader(read):
    """Name of the leader trimmed from this pair before alignment, or None."""
    name = read.query_name
    if not name.startswith(SL_TRIMMED):
        return None
    prefix = name.split(":", 1)[0]
    return prefix.split(".", 1)[1] if "." in prefix else "SL"


def sl_aligned(clip, aligned, sl_sequence, min_overlap, max_aligned=3):
    """Number of spliced-leader 3' bases aligned to the genome, or None if the soft clip does not hold this SL.

    The soft clip must end with (a suffix of) the SL.  The SL's last bases often
    match the genome, because the trans-splice acceptor (...AG, ...CAG) resembles
    the SL's own 3' end (...AACAG); then the clip stops up to `max_aligned` bases
    early and those SL bases start the alignment.
    """
    if not sl_sequence or len(clip) < min_overlap:
        return None
    clip, aligned = clip.upper(), aligned.upper()
    for k in range(max_aligned + 1):
        sl = sl_sequence[:len(sl_sequence) - k]
        n = min(len(clip), len(sl))
        if n >= min_overlap and clip[-n:] == sl[-n:] and aligned[:k] == sl_sequence[len(sl_sequence) - k:]:
            return k
    return None


def sl_status(read, clip, aligned, leaders, min_overlap):
    """(SL bases aligned to the genome, leader name) if READ1 carries a spliced leader, else (None, None).

    A READ1 whose SL was trimmed before alignment starts at the trans-splice site (0 bases aligned).
    Otherwise the leaders are looked for in the soft clip (sl_aligned); the longest match wins.
    """
    name = trimmed_leader(read)
    if name is not None:
        return 0, name
    best = (None, None, -1)
    for leader, seq in leaders:
        k = sl_aligned(clip, aligned, seq, min_overlap)
        if k is not None and len(seq) - k > best[2]:
            best = (k, leader, len(seq) - k)
    return best[0], best[1]


def sl_trimmed_5p(read, pos, strand, clip):
    """5' end of an SL-trimmed READ1, with a short 5' soft clip added back.

    After the SL (and the extra G before it) is trimmed, the read starts at the trans-splice site and holds no
    non-templated bases; first bases that differ from the genome (sequencing errors, strain variants) must not move
    the CTSS downstream when the aligner soft-clipped them.
    """
    if read.query_name.startswith(SL_TRIMMED) and 0 < len(clip) <= MAX_SL_CLIP:
        new = pos - len(clip) if strand == "+" else pos + len(clip)
        if 0 <= new < read.header.get_reference_length(read.reference_name):  # not past a contig end
            return new
    return pos


def sl_trimmed_unplaceable(read, clip):
    """SL-trimmed READ1 with a longer 5' soft clip: its 5' end cannot be placed (extra non-genomic bases after the
    leader, or an indel near the 5' end), so it gives no TSS."""
    return read.query_name.startswith(SL_TRIMMED) and len(clip) > MAX_SL_CLIP


def first_g_mismatch(read):
    """READ1 aligned from its first 5' base, which is a G that differs from the genome (MD tag): the non-templated
    G added by the reverse transcriptase at the cap; the TSS is the next base (as CAGEr's removeFirstG)."""
    if read.query_name.startswith(SL_TRIMMED) or not read.has_tag("MD"):
        return False
    cigar = [op for op in (read.cigartuples or []) if op[0] != 5]
    seq = read.query_sequence or ""
    if not cigar or not seq:
        return False
    md = read.get_tag("MD")
    if read.is_reverse:
        return cigar[-1][0] in (0, 7, 8) and seq[-1] == "C" and bool(MD_LAST.search(md))
    return cigar[0][0] in (0, 7, 8) and seq[0] == "G" and bool(MD_FIRST.match(md))
