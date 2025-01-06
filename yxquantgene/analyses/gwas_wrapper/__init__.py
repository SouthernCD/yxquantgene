from yxquantgene.analyses.gwas_wrapper.magma_wrapper import MAGMA_JOB
from yxquantgene.analyses.gwas_wrapper.gemma_wrapper import GEMMA_JOB

if __name__ == '__main__':
    import pandas as pd
    from yxutil import mkdir

    # input files
    
    # input vcf file
    # A VCF file with genotypes, sample order should be the same as other files
    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/ai.vcf"
    # phenotype_file
    # A single column CSV file without header, values as phenotype, sample order should be the same as vcf_file
    phenotype_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/phenotype.csv"
    # kinship_file
    # A CSV matrix without header and index, values as kinship values, sample order should be the same as phenotype_file
    kinship_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/kinship.csv"
    # gff_file
    # A standard GFF3 file
    gff_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/Sbicolor_313_v3.1.gene.gff3"

    # running
    # GEMMA_JOB
    name = "ai"
    work_dir = "ai_gemma"
    mkdir(work_dir)  

    gemma_job = GEMMA_Job(name=name, work_dir=work_dir, phenotype_file=phenotype_file, input_vcf_file=input_vcf_file, kinship_file=kinship_file)
    gemma_job.run()
    gemma_results_df = pd.read_hdf(gemma_job.gwas_df_h5, 'df')

    # output files
    
    # output gwas df h5 file
    # A HDF5 file with key 'df' and columns: 
    ## For site gwas, gwas_df should have columns: ID, CHROM, POS, REF, ALT, BETA, P_VALUE, FDR
    gemma_results_df

    # MAGMA_JOB
    magma_job = MAGMA_JOB(name=name, work_dir=work_dir, input_vcf_file=input_vcf_file,
                          gff_file=gff_file, site_gwas_df_h5=gemma_job.gwas_df_h5)
    magma_job.run()
    magma_results_df = pd.read_hdf(magma_job.gwas_df_h5, 'df')

    # output files
    # A HDF5 file with key 'df' and columns:
    ## For gene gwas, gwas_df should have columns: ID, CHROM, START, END, NSNPS, P_VALUE, FDR


