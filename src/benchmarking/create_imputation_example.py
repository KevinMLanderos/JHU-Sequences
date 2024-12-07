import sys
from pathlib import Path
import pandas as pd
import numpy as np

"""
Roughly follow Beagle's protocol for generating benchmarking examples with some extra code for ensuring we have ground truth.
Note this will take GBs of RAM for large reference panels. This was run on a compute node on Rockfish.
Author: Nicole Brown
"""

# HELPERS
def save_vcf(header_lines, vcf_info_df, output_file):
    # remove output_file if it exists
    output_file = Path(output_file)
    if output_file.exists():
        output_file.unlink()

    # Write header
    with open(output_file, 'w') as f:
        for comment in header_lines:
            f.write(comment + "\n")

    # Append df to the file
    vcf_info_df.to_csv(output_file, index=False, sep='\t', mode='a')
    return

# MAIN
def create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder):
    output_folder = Path(output_folder)

    # Read only the header lines (lines starting with '##' or '#')
    with open(input_vcf_file, 'r') as file:
        header_lines = []
        column_line = None
        for line in file:
            if line.startswith('##'): # random info at the top
                header_lines.append(line.strip())
            elif line.startswith('#'):
                column_line = line.strip() # Get column names
                break

    column_names = column_line.split("\t")
    info_column_names = column_names[0:9]
    sample_column_names = column_names[9:max_samples+9]
    vcf_info = pd.read_csv(input_vcf_file, sep='\t', comment='#', header=None, names=column_names, nrows=max_positions, usecols=column_names[0:max_samples+9])

    # convert to phased, GT Format - our code cannot use likelihoods so we want to compare apples to apples
    # also impute expects phased reference panels
    vcf_info["FORMAT"] = "GT"
    vcf_info[sample_column_names] = vcf_info[sample_column_names].apply(lambda x: x.str[0] + "|" + x.str[2])

    # drop duplicate rows - beagle expects one record per position and crashes otherwise
    vcf_info = vcf_info.drop_duplicates(subset='POS', keep='first')

    split_index = int(len(sample_column_names) * 0.6)
    samples_ref = info_column_names + sample_column_names[:split_index]
    samples_target = info_column_names + sample_column_names[split_index:]

    # take 60% of these columns - these are the reference cols now
    save_vcf(header_lines, vcf_info[samples_ref], output_folder / f"ref_{max_samples}_{max_positions}.vcf")
    # other 30% are target cols
    save_vcf(header_lines, vcf_info[samples_target], output_folder / f"target_gt_{max_samples}_{max_positions}.vcf")

    # tell user total rows and columns in target
    # tell user number of multiallelic sites
    target = vcf_info[samples_target]
    print("Target shape", target.shape, flush=True)
    multiallelic_sites = target.loc[target['ALT'].str.contains(',')]
    print("# Multi-allelic sites:", multiallelic_sites.shape[0], flush=True)

    # remove all haplotypes from multiallelic sites - this is primarily what we want to test
    num_samples = target.shape[1] - 9
    total_sites_to_scramble = target.shape[0] * num_samples * .3
    total_multiallelic_sites = multiallelic_sites.shape[0] * num_samples
    total_to_randomly_choose = int(total_sites_to_scramble - total_multiallelic_sites)
    target.loc[multiallelic_sites.index, target.columns[9:]] = "./."

    # remove a couple random other sites - arbitrarily, lets create 30% missing spots to impute
    random_rows = np.random.choice(target.index, total_to_randomly_choose, replace=True).astype(int)
    random_cols = np.random.choice(target.columns[9:], total_to_randomly_choose, replace=True)

    for i in range(len(random_rows)):
        target.loc[random_rows[i], random_cols[i]] = "./."

    # save target
    save_vcf(header_lines, target, output_folder / f"target_{max_samples}_{max_positions}.vcf")
    return

# toy example derived from beagle
input_vcf_file = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/beagle/test.vcf"
output_folder = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/beagle/"
max_samples = 100
max_positions = 1000
create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder)

# chr 14 - 1,000 positions, 1000 samples
input_vcf_file = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/chr14.vcf"
output_folder = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"
max_samples = 1000
max_positions = 1000
create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder)

# chr 14 - 10,000 positions, 1000 samples
input_vcf_file = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/chr14.vcf"
output_folder = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"
max_samples = 1000
max_positions = 10000
create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder)

# chr 14 - 100,000 positions, 1000 samples
input_vcf_file = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/chr14.vcf"
output_folder = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"
max_samples = 1000
max_positions = 100000
create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder)

# chr 14 - 1M, 1000 samples
input_vcf_file = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/chr14.vcf"
output_folder = "/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"
max_samples = 1000
max_positions = 1000000
create_imputation_example(input_vcf_file, max_samples, max_positions, output_folder)
