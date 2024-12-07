# Gene Imputation for Multi-Allelic Sites
Applications of the Positional Burrow Wheeler Transform (PBWT)   
Members: Yi Yang, Kevin Meza Landeros, Kiki Zhang, Nicole Brown

-----

Our method:  
To run our algorithm, input files should be obtained from the data folder (reference VCF: data/ref.vcf, target VCF: data/target.vcf), which constitute a small subset of variants from chromosome 22. Also python dependencies should be installed, particularly pysam, which is used to handle VCF files. Script “src/imputation_with_PBWT 2.py” should be run, and it will produce a VCF file as an output containing the imputed alleles named “imputed_results_updated.vcf”.

Parallelization evaluation analysis: run src/parallelization_with_mPBWT, it would run the parallelization, evaluate the mPBWT in both parallel and sequential across different scenarios, and output the comparison plot. 

Benchmarking Analysis: All code for benchmarking can be found in the Variant_calling_PBWT/src/benchmarking/ folder.

Benchmarking first requires downloading chromosome 14 from the 1000 genomes project, phase 3. This can be found here: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ (ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz).
Once this file is downloaded, you can run create_imputation_example.py and create_imputation_example_impute.py with python to generate the examples for all tools to run on. You will need to adjust the output folder file paths starting on line 90 of create_imputation_example.py and line 112 of create_imputation_example_impute.py to your chosen file paths (any file path will work). You should see multiple output folders with  *target_gt.vcf, *target.vcf, and *ref.vcf files generated.
Once you have the examples generated, simply run run_benchmark_beagle.sh, run_benchmark_glimpse.sh, and run_benchmark_impute.sh to run benchmarks for each tool. The *_benchmark.txt files that are generated will contain benchmarking information. Again, adjust the file_path_prefix to point to the parent folder where  *target_gt.vcf, *target.vcf, and *ref.vcf are found, and set the threads parameter according to how many threads you would like to run. Make sure your output directory exists before running the script. The output directory should contain a *imputed.vcf file with imputation results for the tool tested.
Post-run, you can compute accuracy for GLIMPSE and Beagle results using get_imputation_accuracy.py (accuracy at all missing sites) and get_imputation_accuracy_ma.py (accuracy at just multi-allelic sites). Both scripts take the *target_gt.vcf, *target.vcf, and *imputed.vcf files as positional arguments 1, 2, and 3 on the command line.

Credits: 
Our method is inspired and based on the PBWT originally published by Durbin et al. in 2014 and mPBWT published by Ardlan et al. in 2019
Durbin, R. (2014). Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT). Bioinformatics, 30(9), 1266-1272.
Naseri, A., Zhi, D., & Zhang, S. (2019). Multi-allelic positional Burrows-Wheeler transform. BMC bioinformatics, 20, 1-8.
