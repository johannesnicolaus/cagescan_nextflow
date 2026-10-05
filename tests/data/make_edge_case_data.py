#!/usr/bin/env python3
# Edge cases for StringTie on CAGEscan data (source of tests/data/modules/edge_*), reusing make_test_data.py.
# usage: python3 make_edge_case_data.py OUTDIR
import importlib.util
import os
import sys

spec = importlib.util.spec_from_file_location("sim", os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_test_data.py"))
sim = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim)
sim.GENES = [
    # A: two TSSs inside the same first exon (200 nt apart), same downstream exons
    ("chr1", "geneA", "+", {
        "geneA.up":   ([(5001, 5600), (6001, 6400), (7001, 8200)], 0.5),
        "geneA.down": ([(5201, 5600), (6001, 6400), (7001, 8200)], 0.5)}, False),
    # B: exon skipping close to the TSS (inside insert size)
    ("chr1", "geneB", "+", {
        "geneB.1": ([(20001, 20200), (21001, 21100), (22001, 23500)], 0.6),
        "geneB.2": ([(20001, 20200), (22001, 23500)], 0.4)}, False),
    # C: exon skipping far from the TSS (~1.3 kb downstream, beyond most fragments)
    ("chr1", "geneC", "-", {
        "geneC.1": ([(40001, 41000), (42001, 42150), (43001, 44200), (45001, 45200)], 0.5),
        "geneC.2": ([(40001, 41000), (43001, 44200), (45001, 45200)], 0.5)}, False),
    # D: two independent alternative events (skip exon 2, skip exon 4) -> 4 isoforms, unequal
    ("chr2", "geneD", "+", {
        "geneD.11": ([(5001, 5150), (6001, 6080), (7001, 7120), (8001, 8070), (9001, 10500)], 0.4),
        "geneD.01": ([(5001, 5150), (7001, 7120), (8001, 8070), (9001, 10500)], 0.1),
        "geneD.10": ([(5001, 5150), (6001, 6080), (7001, 7120), (9001, 10500)], 0.1),
        "geneD.00": ([(5001, 5150), (7001, 7120), (9001, 10500)], 0.4)}, False),
]
sim.CHROM_LEN = {"chr1": 50000, "chr2": 15000}
sys.argv = ["x", "--outdir", sys.argv[1], "--pairs", "6000", "--seed", "7"]
sim.main()
