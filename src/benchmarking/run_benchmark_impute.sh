#!/bin/bash
# Loop through increasingly large benckmark files and time

set -e

refs=("ref_1000_1000.vcf" "ref_1000_10000.vcf" "ref_1000_100000.vcf" "ref_1000_1000000.vcf")
gts=("target_1000_1000.vcf" "target_1000_10000.vcf" "target_1000_100000.vcf" "target_1000_1000000.vcf")
outs=("results/impute_10/imputed_1000_1000" "results/impute_10/imputed_1000_10000" "results/impute_10/imputed_1000_100000" "results/impute_10/imputed_1000_1000000")
file_path_prefix="/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes_impute/"

for i in "${!refs[@]}"; do
    bgzip -c "${file_path_prefix}${refs[i]}" > "${file_path_prefix}${refs[i]}.gz"
    bgzip -c "${file_path_prefix}${gts[i]}" > "${file_path_prefix}${gts[i]}.gz"
    bcftools index "${file_path_prefix}${refs[i]}.gz"
    bcftools index "${file_path_prefix}${gts[i]}.gz"
    # placeholder instead of ./.; impute doesn't handle missing alleles
    # bcftools +setGT "${file_path_prefix}${gts[i]}.gz" -Oz -o "${file_path_prefix}${gts[i]}_fixed.gz" -- -t '.' -n '0'
    # bcftools index "${file_path_prefix}${gts[i]}_fixed.gz"

    /usr/bin/time -v /home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/impute/impute5_v1.2.0/impute5_v1.2.0_static \
     --no-out-gp-field --no-out-ds-field --no-out-index \
     --h "${file_path_prefix}${refs[i]}.gz" \
     --g "${file_path_prefix}${gts[i]}.gz" \
     --o "${file_path_prefix}${outs[i]}.vcf.gz" \
     --r 14:0-52670879 \
     --buffer-region 14:0-52670879 > "${file_path_prefix}${outs[i]}_mylog.txt" 2> "${file_path_prefix}${outs[i]}_benchmark.txt" \
     --threads 10
    #  --haploid

    gunzip -f "${file_path_prefix}${outs[i]}.vcf.gz"
    echo "Impute run completed for ref=${refs[i]}, gt=${gts[i]}, out=${outs[i]}.vcf"
done

echo "All Impute runs completed!"
