#!/bin/bash
# Loop through increasingly large benckmark files and time

set -e

refs=("ref_1000_1000.vcf" "ref_1000_10000.vcf" "ref_1000_100000.vcf" "ref_1000_1000000.vcf")
gts=("target_1000_1000.vcf" "target_1000_10000.vcf" "target_1000_100000.vcf" "target_1000_1000000.vcf")
outs=("results/beagle_10/imputed_1000_1000" "results/beagle_10/imputed_1000_10000" "results/beagle_10/imputed_1000_100000" "results/beagle_10/imputed_1000_1000000")
file_path_prefix="/home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/data/1000genomes/"

for i in "${!refs[@]}"; do
    /usr/bin/time -v java -jar /home/nbrown62/data_mschatz1/nbrown62/seqs/solutions/final_project/beagle/beagle.jar \
        ref="${file_path_prefix}${refs[i]}" \
        gt="${file_path_prefix}${gts[i]}" \
        nthreads=10 \
        out="${file_path_prefix}${outs[i]}" >  "${file_path_prefix}${outs[i]}_mylog.txt" 2> "${file_path_prefix}${outs[i]}_benchmark.txt"
    gunzip -f "${file_path_prefix}${outs[i]}.vcf.gz"
    echo "Beagle run completed for ref=${refs[i]}, gt=${gts[i]}, out=${outs[i]}"
done

echo "All Beagle runs completed!"
