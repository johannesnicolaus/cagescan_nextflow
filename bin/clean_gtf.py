#!/usr/bin/env python3
"""Drop GTF records with an empty transcript_id.

NCBI-style GTFs (e.g. RefSeq / GenBank annotation exports) give gene rows `transcript_id ""`; gffcompare (and gffread)
reject such records ("no valid ID found for GFF record"). Gene rows are not needed downstream (STAR, StringTie, the
pipeline's scripts and gffcompare use transcript and exon records), so they are dropped; all other lines are copied
unchanged, so annotations without empty transcript IDs are unchanged.

usage: clean_gtf.py in.gtf out.gtf
"""

import re
import sys

EMPTY_TRANSCRIPT_ID = re.compile(r'(^|;)\s*transcript_id\s+""')

src, dst = sys.argv[1], sys.argv[2]
dropped = 0
with open(src) as handle, open(dst, "w") as out:
    for line in handle:
        fields = line.split("\t")
        if not line.startswith("#") and len(fields) >= 9 and EMPTY_TRANSCRIPT_ID.search(fields[8]):
            dropped += 1
            continue
        out.write(line)
print(f"{src}: dropped {dropped} record(s) with an empty transcript_id", file=sys.stderr)
