# Credits: 
# Our method is inspired and based on the PBWT originally published by Durbin et al. in 2014 and mPBWT published by Ardlan et al. in 2019
# Durbin, R. (2014). Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT). Bioinformatics, 30(9), 1266-1272.
# Naseri, A., Zhi, D., & Zhang, S. (2019). Multi-allelic positional Burrows-Wheeler transform. BMC bioinformatics, 20, 1-8.


from collections import Counter
import random
import pysam
import pandas as pd
import numpy as np

# Convert VCFs to matrixes
def vcf_to_haplotype_matrices(vcf_file):
    """
    Function that reads a VCF file and produces two matrices:
    one for maternal alleles and one for paternal alleles.
    Rows are samples, columns are variants, and values represent:
    0 for the reference allele, 1, 2, ... for alternate alleles, and NaN for missing data.
    """
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)
    # Extract positions and samples
    positions = []
    samples = list(vcf.header.samples)
    # Initialize dictionaries to store matrix data for maternal and paternal alleles
    maternal_data = {sample: [] for sample in samples}
    paternal_data = {sample: [] for sample in samples}
    for record in vcf.fetch():
        positions.append(record.pos)  # Save the genomic position
        # Process each sample
        for sample in samples:
            genotype = record.samples[sample]["GT"]
            # Handle missing data
            if genotype is None:
                maternal_data[sample].append(-1)
                paternal_data[sample].append(-1)
            else:
                # Assign maternal and paternal alleles
                maternal, paternal = genotype
                maternal_value = -1 if maternal is None else maternal # (maternal if maternal != 0 else 0)
                paternal_value = -1 if paternal is None else paternal # (maternal if maternal != 0 else 0)
                # Append values to respective matrices
                maternal_data[sample].append(maternal_value)
                paternal_data[sample].append(paternal_value)
    # Convert dictionaries to DataFrames
    maternal_df = pd.DataFrame(maternal_data, index=positions).transpose()
    paternal_df = pd.DataFrame(paternal_data, index=positions).transpose()
    # Return the two matrices
    return maternal_df, paternal_df

# Algorithm 3: find long matches
def find_long_matches(xk, t, ak, dk, k, L):
    """
    Reports haplotype matches with start and end positions.
    Parameters:
        xk (list): Contains alleles at position k for all haplotypes.
        t (int): Number of possible alleles (e.g., t-allelic sites).
        ak (list): PBWT ordering of indices at position k.
        dk (list): Divergence array at position k.
        k (int): Current position index.
        L (int): Minimum match length to consider.
    Outputs:
        Reports matches between haplotypes, including start and end positions.
    """
    matches = []
    # Initialize variables
    m = [False] * t  # Tracks the presence of alleles
    i0 = 0           # Start index of the current interval
    M = len(dk)      # Total number of haplotypes
    # Iterate through the haplotypes in PBWT order
    for i in range(M):
        # Check if divergence exceeds the threshold (report when # of missmatches (divergence) is big)
        if dk[i] > k - L: # NOTE: Shall I add equal
            # Identify distinct alleles in the current interval
            distinct_alleles = [j for j in range(t) if m[j]]
            if len(distinct_alleles) >= 2:  # At least two distinct alleles
                # Process matches in the interval
                for ia in range(i0, i):
                    dmin = 0
                    for ib in range(ia + 1, i):
                        if dk[ib] > dmin:
                            dmin = dk[ib]
                        if xk[ak[ia]] == xk[ak[ib]] and any([ak[ia] == (len(ak)-1), ak[ib] == (len(ak)-1)]):  # If alleles match 
                            print(f"Match: Haplotype {ak[ia]} and Haplotype {ak[ib]}, Start: {dmin}, End: {k}, Len: {k-dmin}")
                            matches.append((ak[ia], ak[ib], dmin, k, k-dmin))
                        else: # xk[ak[ia]] != xk[ak[ib]]
                            if k - dmin >= L and any([ak[ia] == (len(ak)-1), ak[ib] == (len(ak)-1)]):  # Ensure the match length meets the threshold 
                                print(f"Match: Haplotype {ak[ia]} and Haplotype {ak[ib]}, Start: {dmin}, End: {k}, Len: {k-dmin}")
                                matches.append((ak[ia], ak[ib], dmin, k, k-dmin))
            # Reset the start index and allele presence tracking
            i0 = i
            m = [False] * t
        # Update allele presence for the current haplotype
        allele = xk[ak[i]]
        m[allele] = True
    # Final check for the remaining interval
    distinct_alleles = [j for j in range(t) if m[j]]
    if len(distinct_alleles) >= 2:  # At least two distinct alleles
        for ia in range(i0, M):
            dmin = 0
            for ib in range(ia + 1, M):
                if dk[ib] > dmin:
                    dmin = dk[ib]
                if xk[ak[ia]] == xk[ak[ib]] and any([ak[ia] == (len(ak)-1), ak[ib] == (len(ak)-1)]):  # If alleles match 
                    print(f"Match: Haplotype {ak[ia]} and {ak[ib]}, Start: {dmin}, End: {k}, Len: {k-dmin}")
                    matches.append((ak[ia], ak[ib], dmin, k, k-dmin))
                else: # xk[ak[ia]] != xk[ak[ib]]
                    if k - dmin >= L and any([ak[ia] == (len(ak)-1), ak[ib] == (len(ak)-1)]):  # Ensure the match length meets the threshold 
                        print(f"Match: Haplotype {ak[ia]} and {ak[ib]}, Start: {dmin}, End: {k}, Len: {k-dmin}")
                        matches.append((ak[ia], ak[ib], dmin, k, k-dmin))
    return matches

# Algorithm 2: Update prefix and divergence arrays
def prefix_and_divergence(xk, t, ak, dk, k): # does not consider an input s
    """
    Implements the prefix_and_divergence procedure.
    Parameters:
        xk (list): List of t-allelic variant sites at position k.
        t (int): Number of alleles.
        ak (list): PBWT ordering of indices at position k.
        dk (list): Divergence array at position k.
        k (int): Current position index.
    Returns:
        ak_next (list): Updated PBWT ordering for position k+1.
        dk_next (list): Updated divergence array for position k+1.
    """
    # Initialize arrays for storing a[t], d[t], p[t], and u[t]
    a = [[] for _ in range(t)]  # a[t] will store arrays for each allele
    d = [[] for _ in range(t)]  # d[t] will store divergence arrays for each allele
    p = [k + 1] * t  # Initialize p[j] to k + 1 for all j
    u = [0] * t      # Initialize u[j] to 0 for all j
    M = len(ak)  # Number of haplotypes
    # Loop through all haplotypes
    for i in range(M):
        allele = xk[ak[i]]  # Identify the allele for the current haplotype
        # Append the current haplotype index to the appropriate a[allele]
        a[allele].append(ak[i])
        # Update p[j] values for all alleles
        for j in range(t):
            if dk[i] > p[j]:
                p[j] = dk[i]
        # Store the updated p[allele] in d[allele]
        d[allele].append(p[allele])
        # Reset p[allele] for the current allele
        p[allele] = 0
        # Increment the counter for the current allele
        u[allele] += 1
    # Concatenate arrays a[i] and d[i] for all alleles to form ak+1 and dk+1
    ak_next = [item for sublist in a for item in sublist]
    dk_next = [item for sublist in d for item in sublist]
    return ak_next, dk_next

def imputation(matches, k, X):
    """"
    Calculates consensus allele from previous matches with the highest lenght. Whenever 2 alleles with the same max Len are equally soported, an allele is chosen randomly;
    for example this is the case when only the previopus allele is matching and there are equal number of 0's and 1's.
    Parameters:
        matches (list): List of lists with structure: [Haplotype1, Haplotype2, start, End, k, Len]
        k (int): Current position index.
    Returns:
        consensus_allele (int): Consensus allele for quey sequence in next position.
    """
    max_value = 0; alleles = []
    # Iterate over the list of lists once
    for match in matches:
        value = match[4] # focus on match lenght
        if value > max_value:
            max_value = value
            print(match[:2], len(X))
            hap = [x for x in match[:2] if x != len(X)-1][0]; print(hap)
            alleles = [X[hap][k]]  # Reset to the new max
        elif value == max_value:
            print(match[:2], len(X))
            hap = [x for x in match[:2] if x != len(X)-1][0]; print(hap)
            alleles.append(X[hap][k])
    # Get all numbers with the maximum count
    counts = Counter(alleles); max_count = max(counts.values())
    print(f"counts: {counts}")
    candidates = [number for number, count in counts.items() if count == max_count]
    print(f"candidates: {candidates}")
    # Gewt consensus. Choose one randomly if there's a tie
    consensus_allele = random.choice(candidates)
    return consensus_allele

def iterate_over_positions(s, X, L, t):
    """
    Combines prefix_and_divergence, find_long_matches
    and imputation to iterate over all positions k.
    Parameters:
        s (list): New sequence to match. WILL ALWAYS BE APPENDED AT THE END OF THE REFERENCE MATRIX
        X (list of lists): Reference panel of haplotypes (list of sequences).
        L (int): Minimum match length.
        t (int): Number of alleles.
    Outputs:
        Matches (list) Imputed alles for input query
    """
    # Include new sequence to  X
    X = X + [s] 
    # Initialize PBWT ordering and divergence arrays
    ak = list(range(len(X)))  # Initial PBWT order
    dk = [0] * len(X)         # Initial divergence array
    M = len(X)                # Number of sequences in the reference panel
    N = len(s)                # Length of the sequences
    # Iterate over positions
    for k in range(N):
        matches_=[]
        print(f"\nPosition {k}")
        # Extract xk for the current position from the reference panel (used to get k+1 order)
        xk_current = [X[m][k] for m in range(len(X))]
        print(f"xk_current: {xk_current}")
        if k!=0: # NOTE: Run only if X[M-1][k] == -1 
            # Extract xk for the previous position from the reference panel (used to get matches)
            xk_prev = [X[hap][k-1] for hap in ak]
            print(f"xk_prev: {xk_prev}")
            # Step 1: Find long matches for the current position
            print("Finding matches...")
            matches_ = find_long_matches(xk_prev, t, ak, dk, k, L)
            print(matches_)
        # Step 2: Impute query's next position, only if missing
        if X[M-1][k] == -1:  
            print("Imputing missing allele...")
            # if k==0: # NOTE: Set missing variant to 0 if missing allele in on first position. Might need adjustment in future, and use following variants to impute
            #     imp_allele=0
            if all(i == 0 for i in dk[1:-1]): # NOTE: Modify in the future for first position
                imp_allele=max(set(xk_current), key=xk_current.count)
            else:
                if matches_:
                    print("Calling Function")
                    imp_allele=imputation(matches_, k, X)
                else: # If list is empty # NOTE: Check in future
                    imp_allele=xk_current[0]
            X[M-1][k]=imp_allele # Assign value in query for next position
            xk_current = [X[m][k] for m in range(len(X))]
            print(f"Imputed allele: {imp_allele}")
            print(f"xk_current (after imputation): {xk_current}")
            if imp_allele == -1:
                return "bad"
        # Step 3: Obtain Prefix and divergence arrays for k+1
        print("Calculating prefix and divergence array for next position...")
        ak, dk = prefix_and_divergence(xk_current, t, ak, dk, k)
        print(ak, dk)
    return X[M-1] # matches_ # X[M-1] 


###### ------------------ Main ------------------ ######

# - Define parameters
t = 2  # Biallelic sites
L = 1  # Minimum match length

# - Set reference
vcf_file="/Users/kikizhang/Downloads/pbwt/ref.vcf"
ref_m1, ref_m2 = vcf_to_haplotype_matrices(vcf_file)
ref_m1 = ref_m1.values.tolist(); ref_m2 = ref_m2.values.tolist() # Convert Df to list of lists

# - Set query
vcf_file="/Users/kikizhang/Downloads/pbwt/target.vcf"
query_m1, query_m2 = vcf_to_haplotype_matrices(vcf_file)
# Check unique values in query
# m1
unique_values = query_m1.stack().unique()
print(unique_values)
# m2
unique_values = query_m2.stack().unique()
print(unique_values)


# - Impute maternal alleles (m1)
print("Imputing maternal alleles...")
imputed_m1 =[]
for index, query in enumerate(query_m1.values.tolist()):
    # imputed_m1 += [iterate_over_positions(s=query, X=ref_m1, L=L, t=t)]
    print(f"# ----- Query #{index} ----- #") 
    r = iterate_over_positions(s=query, X=ref_m1, L=L, t=t)
    if r == "bad":
        break
    else:
        imputed_m1 += [r]

print("DONE.")

# - Impute maternal alleles (m2)
print("Imputing maternal alleles...")
imputed_m2 =[]
for index, query in enumerate(query_m2.values.tolist()): 
    print(f"# ----- Query #{index} ----- #") 
    imputed_m2 += [iterate_over_positions(s=query, X=ref_m1, L=L, t=t)]

print("DONE.")

def imputed_to_vcf(imputed_m1, imputed_m2, original_vcf, output_file):
    """
    Convert imputed sequences back to VCF format and force phasing with text replacement
    """
    # Create VCF as before
    vcf = pysam.VariantFile(original_vcf)
    header = vcf.header
    out_vcf = pysam.VariantFile(output_file + ".temp", 'w', header=header)
    samples = list(vcf.header.samples)
    
    for idx, record in enumerate(vcf.fetch()):
        new_record = out_vcf.new_record(
            contig=record.contig,
            start=record.start,
            stop=record.stop,
            alleles=record.alleles,
            id=record.id,
            qual=record.qual,
            filter=record.filter,
            info=dict(record.info)
        )
        
        for sample_idx, sample in enumerate(samples):
            maternal = imputed_m1[sample_idx][idx]
            paternal = imputed_m2[sample_idx][idx]
            
             # Create new genotype
            new_record.samples[sample]['GT'] = (maternal, paternal)
            
            # Copy any other format fields from original record
            for field in record.format:
                if field != 'GT':
                    new_record.samples[sample][field] = record.samples[sample][field]
        
            
        out_vcf.write(new_record)
    
    out_vcf.close()

    vcf.close()
       # Replace all '/' with '|' in the output file
    with open(output_file + ".temp", 'r') as f:
        content = f.read()
    
    with open(output_file, 'w') as f:
        f.write(content.replace('/', '|'))
    
    # Remove temporary file
    import os
    os.remove(output_file + ".temp")


# Run the VCF conversion
imputed_to_vcf(
    imputed_m1=imputed_m1,
    imputed_m2=imputed_m2,
    original_vcf="/Users/kikizhang/Downloads/pbwt/target.vcf",
    output_file="imputed_results_updated.vcf"
)


######
