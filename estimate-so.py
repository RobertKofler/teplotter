#!/usr/bin/env python
import argparse
import modules


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

args = parser.parse_args()
writer=modules.Writer(args.outfile)

# than normalize each entry
for se in modules.SeqEntryReader(args.seqentry):
    avcov=modules.NormFactor.computeNormFactorForSe([se],args.enddist,args.quantile)
    topr=[se.seqname,f"{avcov:.2f}",str(len(se.cov))]
    tp="\t".join(topr)
    writer.write(tp)
