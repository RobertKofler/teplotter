#!/usr/bin/env python
import argparse
from modules import SeqEntryReader, NormFactor, Writer




parser = argparse.ArgumentParser(description="""           
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--so', type=str, default=None,dest="seqentry", required=True, help="A sequence overview (so) file")
parser.add_argument("--scg-end", type=str, required=False, dest="scgend", default="_scg", help="the ending by which to recognize single copy gens (scg)")
parser.add_argument("--end-distance", type=int, required=False, dest="enddist", default=100, help="distance from ends for normalizing")
parser.add_argument("--exclude-quantile", type=int, required=False, dest="quantile", default=25, help="exclude the most extreme coverage quantiles for normalizing")
parser.add_argument("--output-file", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")

args = parser.parse_args()
writer = Writer(args.outfile)
# first get the normalization factor
normfactor = NormFactor.getNormalizationFactor(args.seqentry, args.scgend ,args.enddist, args.quantile)

# than normalize each entry
for se in SeqEntryReader(args.seqentry):
    sen = se.normalize(normfactor)
    writer.write(str(sen))
