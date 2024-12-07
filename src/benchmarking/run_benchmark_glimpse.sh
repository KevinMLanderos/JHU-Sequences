#!/bin/bash
# Loop through increasingly large benckmark files and time

set -e

refs=("ref_1000_1000.vcf" "ref_1000_10000.vcf" "ref_1000_100000.vcf" "ref_1000_1000000.vcf")
gts=("target_1000_1000.vcf" "target_1000_10000.vcf" "target_1000_100000.vcf" "target_1000_1000000.vcf")
outs=("results/glimpse/imputed_1000_1000" "results/glimpse/imputed_1000_10000" "results/glimpse/imputed_1000_100000" "results/glimpse/imputed_1000_1000000")
file_path_prefix="/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"

for i in "${!refs[@]}"; do
    awk 'BEGIN {OFS="\t"} {if($1 ~ /^#/ || NF < 8) print $0; else $8="."; print $0}' "${file_path_prefix}${refs[i]}" > "${file_path_prefix}${refs[i]}_fixed.vcf"
    awk 'BEGIN {OFS="\t"} {if($1 ~ /^#/ || NF < 8) print $0; else $8="."; print $0}' "${file_path_prefix}${gts[i]}" > "${file_path_prefix}${gts[i]}_fixed.vcf"
    bcftools +fill-tags "${file_path_prefix}${refs[i]}" -- -t AC,AN > "${file_path_prefix}${refs[i]}_fixed.vcf"
    bgzip -c "${file_path_prefix}${refs[i]}_fixed.vcf" > "${file_path_prefix}${refs[i]}.gz"
    bgzip -c "${file_path_prefix}${gts[i]}_fixed.vcf" > "${file_path_prefix}${gts[i]}.gz"
    bcftools index "${file_path_prefix}${refs[i]}.gz"
    bcftools index "${file_path_prefix}${gts[i]}.gz"
    /usr/bin/time -v /home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/glimpse/GLIMPSE_phase_static \
    --reference "${file_path_prefix}${refs[i]}.gz" \
    --input "${file_path_prefix}${gts[i]}.gz" \
    --output "${file_path_prefix}${outs[i]}.vcf" \
    --input-region 14:0-52670879 --output-region 14:0-52670879 \
    --thread 1 \
     --log "${file_path_prefix}${outs[i]}_glimpse.log" >  "${file_path_prefix}${outs[i]}_mylog.txt" 2> "${file_path_prefix}${outs[i]}_benchmark.txt"
    echo "GLIMPSE run completed for ref=${refs[i]}, gt=${gts[i]}, out=${outs[i]}"
done

echo "All GLIMPSE runs completed!"
