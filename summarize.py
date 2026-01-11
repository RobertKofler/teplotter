#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import math
import collections
import fileinput
from modules import fasta_loader
prog = re.compile(r"(\d+)([MISDHN])")


parser = argparse.ArgumentParser(description="""           
summarize coverage for diverse features
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--fasta", type=int, required=True, dest="minmq", default=1, help="the fasta file to which reads were mapped")

args = parser.parse_args()
minmq=args.minmq