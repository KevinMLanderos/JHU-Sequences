import sys
import pandas as pd

def read_vcf(input_vcf_file):
    # Read only the header lines (lines starting with '##' or '#')
    with open(input_vcf_file, 'r') as file:
        # header_lines = []
        column_line = None
        for line in file:
            if line.startswith('##'): # random info at the top
                # header_lines.append(line.strip())
                continue
            elif line.startswith('#'):
                column_line = line.strip() # Get column names
                break

    column_names = column_line.split("\t")
    vcf_info = pd.read_csv(input_vcf_file, sep='\t', comment='#', header=None, names=column_names)
    return vcf_info

def calculate_accuracy(gt_vcf_location, target_vcf_location, imputed_vcf_location):
    gt_vcf = read_vcf(gt_vcf_location)
    target_vcf = read_vcf(target_vcf_location)
    imputed_vcf = read_vcf(imputed_vcf_location)

    target_vcf = target_vcf[target_vcf["POS"].isin(imputed_vcf["POS"])].reset_index(drop=True)
    gt_vcf = gt_vcf[gt_vcf["POS"].isin(imputed_vcf["POS"])].reset_index(drop=True)
    imputed_vcf = imputed_vcf[imputed_vcf.columns[9:]]
    target_vcf = target_vcf[imputed_vcf.columns]
    gt_vcf = gt_vcf[imputed_vcf.columns]
    imputed_vcf = imputed_vcf.apply(lambda x: x.str[0] + "|" + x.str[2])

    # get missing sites
    target_vcf_missing = target_vcf[target_vcf == './.']
    imputed_vcf = imputed_vcf[target_vcf == './.']
    gt_vcf = gt_vcf[target_vcf == './.']

    imputed_vcf_accurate = imputed_vcf[imputed_vcf == gt_vcf]
    accuracy = imputed_vcf_accurate.count().sum() / target_vcf_missing.count().sum()
    print("Accuracy is:", accuracy)
    return accuracy * 100

gt_vcf_location = sys.argv[1]
target_vcf_location = sys.argv[2]
imputed_vcf_location = sys.argv[3]
calculate_accuracy(gt_vcf_location, target_vcf_location, imputed_vcf_location)
