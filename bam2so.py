#!/usr/bin/env python
import argparse
import pysam
from modules import SeqBuilder, Writer, load_fasta



parser = argparse.ArgumentParser(description="""           
summarize coverage for diverse features
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--infile', type=str, dest="infile", required=True, help="Input BAM or SAM file path")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="the fasta file to which reads were mapped")
parser.add_argument("--mapqth", type=int, required=False, dest="mapqth", default=5, help="mapping quality threshold; below ambiguous")
parser.add_argument("--mc-snp", type=int, required=False, dest="mcsnp", default=5, help="minimum count of SNPs")
parser.add_argument("--mf-snp", type=float, required=False, dest="mfsnp", default=0.1, help="minimum frequency of SNPs")
parser.add_argument("--mc-indel", type=int, required=False, dest="mcindel", default=3, help="minimum count of indels")
parser.add_argument("--mf-indel", type=float, required=False, dest="mfindel", default=0.01, help="minimum frequency of indels")
parser.add_argument("--output-file", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")

args = parser.parse_args()
writer = Writer(args.outfile)

# load fasta from file into dict
reference_dict = load_fasta(args.fasta)

builder=None

infile_path = args.infile
mode = 'rb' if infile_path.lower().endswith('.bam') else 'r'
samfile = pysam.AlignmentFile(infile_path, mode)

for read in samfile:
    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    if read.is_supplementary:
        continue

    ref_name = read.reference_name
    pos = read.reference_start + 1
    mapq = read.mapping_quality if read.mapping_quality is not None else 0
    cigar = read.cigarstring
    read_sequence = read.query_sequence.upper() if read.query_sequence is not None else ''
    
    if cigar is None or cigar == '*':
        continue
    
    if builder is None:
        ref_sequence = reference_dict[ref_name]
        builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)

    if ref_name != builder.seqname:
        seq_entry = builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
        writer.write(str(seq_entry))
        
        ref_sequence = reference_dict[ref_name]
        builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)

    builder.add_read(pos, cigar, mapq, read_sequence)

samfile.close()

# process the last one as well
seq_entry = None
if builder is not None:
    seq_entry = builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
    writer.write(str(seq_entry))



 