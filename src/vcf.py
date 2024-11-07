import pysam
import pandas as pd
import numpy as np

def vcf_to_haplotype_matrix(vcf_file):
    ''' Function that reads a VCF file, and produce a matrix
    where rows are haplotypes, columns are variants and values
    represent wether the allele is the same ase the reference (0)
    or an alternaive one (1,2,..)'''
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)
    # Extract positions and samples
    positions = []
    samples = list(vcf.header.samples)
    # Initialize a dictionary to store matrix data
    # The matrix dictionary will store columns as positions and rows as haplotypes
    matrix_data = {sample: [] for sample in samples}
    for record in vcf.fetch():
        positions.append(record.pos)  # Save the genomic position
        # Process each sample (haplotype)
        for sample in samples:
            genotype = record.samples[sample]["GT"]
            # Handle missing data
            if genotype is None:
                matrix_data[sample].append(np.nan)
            else:
                # Convert genotypes to matrix values
                # Reference allele is represented by 0
                # Alternative alleles are represented by 1, 2, 3, ...
                if all(allele == 0 for allele in genotype):
                    matrix_data[sample].append(0)  # Matches reference
                else:
                    # Convert to 1-based allele index (1, 2, ...) based on the alt allele
                    allele_number = genotype[0] if genotype[0] != 0 else genotype[1]
                    matrix_data[sample].append(allele_number)
    # Convert to DataFrame with positions as columns and samples as rows
    matrix_df = pd.DataFrame(matrix_data, index=positions).transpose()
    # Return the matrix DataFrame
    return matrix_df

def compute_next_prefix_array(xk, t, ak):
    ''' Function that takes a prefix array (ak) for k-th variant and
    computes the prefix array for the the position (ak+1).'''
    # Step 1: Initialize counters for each allele
    u = [0] * t  # Initialize a counter array with t elements set to 0
    a = [[] for _ in range(t)]  # Initialize an array of lists to store haplotype indices
    # Step 2: Populate the a[allele] arrays based on the current prefix array ak
for i in range(len(ak)):  # i=0 # Iterate through each element in ak
    print(i)
    haplotype_index = ak[i]  # Get the current haplotype index from ak
    allele = xk[haplotype_index]  # Get the allele for this haplotype at position k
    a[allele].append(haplotype_index)  # Append the haplotype index to the respective allele list
    u[allele] += 1  # Increment the counter for this allele
    # Step 3: Concatenate all a[allele] lists to form the next prefix array ak+1
    ak_next = [index for sublist in a for index in sublist]
    return ak_next


############# Main #############

# Generate haplotype-variants matrix
vcf_file = "/Users/kmlanderos/Documents/Johns_Hopkins/Fall_2024/Computational_Genomics_Sequences/Project/data/chr14_SNPs_MULTI_ALLELIC.vcf"  # Replace with the path to your VCF file
haplotype_matrix = vcf_to_haplotype_matrix(vcf_file)
print(haplotype_matrix) # Display the matrix

# xk = haplotype_matrix.iloc[:, 0] # Current column variants
# t = len(np.unique(xk.values)) # # of current column unique variants
# ak = haplotype_matrix.index # Initial prefix array (haplotype order)
# ak_next = compute_next_prefix_array(xk, t, ak)


xk = haplotype_matrix.iloc[:, 114] # Current column variants
t = len(np.unique(xk.values)) # # of current column unique variants
ak = haplotype_matrix.index # Initial prefix array (haplotype order)
ak_next = compute_next_prefix_array(xk, t, ak)


for pos in range(0, haplotype_matrix.shape[0]): # Iterate over every position (columns) in haplotype_matrix
    print(pos)
    xk = haplotype_matrix.iloc[:, pos] # Current column variants
    t = len(np.unique(xk.values)) # # of current column unique variants
    if pos == 0:
        ak = haplotype_matrix.index
    # Get prefix array for next position (ak+1)
    ak_next = compute_next_prefix_array(xk, t, ak)
    ak = ak_next

print(ak)
