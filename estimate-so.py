#!/usr/bin/env python3
import argparse
from modules import SeqEntryReader, Writer, NormFactor


parser = argparse.ArgumentParser(description="""           
estimates the average coverage for each sequence overview entry;
notably, the average coverage corresponds to the copy number if the coverage is normalized to the coverage of single copy genes
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--so', type=str, default=None,dest="seqentry", required=True, help="A sequence overview (so) file")
parser.add_argument("--end-distance", type=int, required=False, dest="enddist", default=100, help="distance from ends for normalizing")
parser.add_argument("--exclude-quantile", type=int, required=False, dest="quantile", default=25, help="exclude the most extreme coverage quantiles for normalizing")
parser.add_argument("--output-file", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")
parser.add_argument("--sample-id", type=str, required=False, dest="sampleid", default="x", help="the ID of current sample")


args = parser.parse_args()
writer = Writer(args.outfile)

# than normalize each entry
for se in SeqEntryReader(args.seqentry):
    cost=NormFactor.getCovStat(se,args.enddist,args.quantile)
    topr=[args.sampleid,se.seqname,str(len(se.cov))]
    for c in cost:
        if c is not None:
            topr.append(f"{c:.2f}",)
        else:
            topr.append("na")
    
    tp="\t".join(topr)
    writer.write(tp)
