#!/programs/bin/R/R

######### R parallelized wrapper for FastEPRR #######################

# You need: vcf files imputed and phased per chromosome
# This wrapper takes 3 arguments:
# 1) vcf_files_path         --> path to the folder with the vcf files, vcf file should be split by chromosome
# 2) window_size            --> rho will be calculated in windows with this size
# 3) chromosome_list_file   --> file with the list of chromosomes to be analyzed
#
# Usage example: Rscript FASTEPRR.R /path/to/vcf_files 100 chromosome_list.txt
#
# Script optimized to work with 40 cores to change this modify line 23, 45 and 46
# Three folders will be created: Step1_output, Step2_output, Step3_output
#
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

library(FastEPRR)
library(foreach)
library(doMC)
registerDoMC(80)

### STEP 1

vcf_files_path <- args[1]
window_size <- as.integer(args[2])
chromosome_list_file <- args[3]
step_size <- window_size/2


# Read chromosome names from the file
chromosome_names <- readLines(chromosome_list_file)

# Create the output folder for Step 1
system("mkdir -p Step1_output")

foreach(chr_name = chromosome_names) %dopar% {
    input <- sprintf("%s/%s.vcf.gz", vcf_files_path, chr_name)
    step1_output <- sprintf("Step1_output/%s.coor", chr_name)
    # # 调试输出
    # cat("Processing chromosome:", chr_name, "\n")
    FastEPRR_VCF_step1(vcfFilePath = input, erStart = 1, winLength = window_size, stepLength = step_size, srcOutputFilePath = step1_output)
}

### STEP 2
system("mkdir -p Step2_output")

foreach(i = 1:80) %dopar% {
    FastEPRR_VCF_step2(srcFolderPath = "Step1_output", jobNumber = 80, currJob = i, DXOutputFolderPath = "Step2_output")
}

### STEP 3
system("mkdir -p Step3_output")
FastEPRR_VCF_step3(srcFolderPath = "Step1_output", DXFolderPath = "Step2_output", finalOutputFolderPath = "Step3_output")

