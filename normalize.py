#!/usr/bin/env python
import argparse
import modules




parser = argparse.ArgumentParser(description="""           
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--seq-entries', type=str, default=None,dest="seqentry", required=True, help="A seq entry file")
parser.add_argument("--scg-end", type=str, required=False, dest="scgend", default="_scg", help="the ending by which to recognize single copy gens (scg)")
parser.add_argument("--end-distance", type=int, required=False, dest="enddist", default=100, help="distance from ends for normalizing")


args = parser.parse_args()
scgend=args.scgend
scgs=[]
for se in modules.SeqEntryReader(args.seqentry):
    print(se.seqname)
    if se.seqname.endswith(scgend):
        scgs.append(se)