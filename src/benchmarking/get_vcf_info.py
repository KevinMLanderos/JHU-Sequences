"""
Tiny utility to check VCF info.

Author: Nicole Brown
"""

import sys
import pandas as pd

input_vcf_file = sys.argv[1]

# Read only the header lines (lines starting with '##' or '#')
with open(input_vcf_file, 'r') as file:
    header_line = None
    num_lines = 0
    for line in file:
        if line.startswith('##'): # random info at the top
            continue
        elif line.startswith('#'):
            header_line = line.strip() # Get column names
        num_lines += 1

print(header_line.split("\t")[0:15])
# Rough sample count - typically there are 9 rows before samples start
print("Num Samples:", len(header_line.split("\t")) - 9)
print("Num Rows:", num_lines)
