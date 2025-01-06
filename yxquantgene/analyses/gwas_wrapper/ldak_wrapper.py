from yxutil import mkdir
import pandas as pd
from yxseq import read_gff_file
from yxutil import cmd_run, have_file
from yxquantgene.utils.vcf import get_chr_id_list_from_vcf, get_sample_list_from_vcf
from cyvcf2 import VCF
import os

LDAK_ENV = 'ldak_env'
LDAK_PATH = 'ldak6'


def get_gentype_file(input_vcf_file, output_genetype_stem):
    """
    ldak expects the genetype file to be in plink format, and chr_id should be in the format of integer, snp id (predictor) should be in the format of chr_id:pos
    """
    chr_id_list = get_chr_id_list_from_vcf(input_vcf_file)
    chr_id_map_dict = {chr_id: i + 1 for i, chr_id in enumerate(chr_id_list)}

    chr_id_map_file = f"{output_genetype_stem}.chr_id_map"
    with open(chr_id_map_file, 'w') as f:
        for old_name, new_name in chr_id_map_dict.items():
            f.write(f"{old_name},{new_name}\n")

    tmp_vcf_file = input_vcf_file.replace(".vcf.gz", ".tmp.vcf")

    vcf_reader = VCF(input_vcf_file)

    # 修改 VCF 文件头中的 contig 命名
    header_lines = vcf_reader.raw_header.split('\n')
    new_header_lines = []
    for line in header_lines:
        if line.startswith('##contig'):
            for old_name, new_name in chr_id_map_dict.items():
                line = line.replace(f'ID={old_name}', f'ID={new_name}')
        new_header_lines.append(line)
    new_header = '\n'.join(new_header_lines)

    with open(tmp_vcf_file, 'w') as f:
        f.write(new_header)

        num = 0
        for record in vcf_reader:
            num += 1
            if num % 1000 == 0:
                print(f"Processed {num} records")

            het_num = record.num_het
            hom_ref_num = record.num_hom_ref
            hom_alt_num = record.num_hom_alt
            maf = (het_num + 2*hom_alt_num)/(2*(hom_ref_num + het_num + hom_alt_num))
            maf = min(maf, 1-maf)
            het_freq = het_num / (het_num + hom_ref_num + hom_alt_num)
            het_exp_freq = 2 * maf * (1 - maf)

            if het_freq > 0.1 or maf < 0.01:
                continue

            record.CHROM = str(chr_id_map_dict[record.CHROM])
            record.ID = f"{record.CHROM}:{record.POS}"

            f.write(str(record))

    vcf_reader.close()

    cmd_run(f'bgzip {tmp_vcf_file}')
    cmd_run(f'tabix -f -p vcf {tmp_vcf_file}.gz')

    cmd = f"plink --vcf {tmp_vcf_file}.gz --make-bed --out {output_genetype_stem}"
    cmd_run(cmd, cwd=os.path.dirname(output_genetype_stem))

    cmd_run(f'rm {tmp_vcf_file}.gz')
    cmd_run(f'rm {tmp_vcf_file}.gz.tbi')

    return output_genetype_stem


def get_summary_statistics_file(gwas_df, sample_number, sum_stat_txt_file, extract_file, chr_id_map_file=None):
    if chr_id_map_file:
        chr_id_map_dict = {}
        with open(chr_id_map_file) as f:
            for line in f:
                old_name, new_name = line.strip().split(',')
                chr_id_map_dict[old_name] = new_name

        gwas_df['CHROM'] = gwas_df['CHROM'].map(chr_id_map_dict)

    sum_stat_df = pd.DataFrame()
    sum_stat_df['Predictor'] = gwas_df[['CHROM', 'POS']].agg(
        lambda x: f"{x['CHROM']}:{x['POS']}", axis=1)
    
    is_snp = (gwas_df['REF'].str.len() == 1) & (gwas_df['ALT'].str.len() == 1)

    sum_stat_df['A1'] = gwas_df['REF'].where(is_snp, 'X')
    sum_stat_df['A2'] = gwas_df['ALT'].where(is_snp, 'Y')
    sum_stat_df['n'] = sample_number
    sum_stat_df['Direction'] = gwas_df['BETA']
    sum_stat_df['P'] = gwas_df['P_VALUE']
    sum_stat_df.to_csv(sum_stat_txt_file, index=False, sep=" ")

    sum_stat_df['Predictor'].to_csv(extract_file, index=False, header=False)

    return sum_stat_txt_file, extract_file


def get_genefile(input_gff3_file, genefile, chr_id_map_file=None):
    if chr_id_map_file:
        chr_id_map_dict = {}
        with open(chr_id_map_file) as f:
            for line in f:
                old_name, new_name = line.strip().split(',')
                chr_id_map_dict[old_name] = new_name

    gene_dict = read_gff_file(input_gff3_file)['gene']
    gene_id_list = []
    for chr_id in chr_id_map_dict:
        chr_gen_list = [
            gene_id for gene_id in gene_dict if gene_dict[gene_id].chr_id == chr_id]
        chr_gen_list = sorted(chr_gen_list, key=lambda x: gene_dict[x].start)
        gene_id_list.extend(chr_gen_list)

    with open(genefile, 'w') as f:
        for gene_id in gene_id_list:
            gene = gene_dict[gene_id]
            if chr_id_map_file:
                if gene.chr_id not in chr_id_map_dict:
                    continue
                f.write(
                    f"{gene_id} {chr_id_map_dict[gene.chr_id]} {gene.start - 1} {gene.end - 1} {gene.strand}\n")
            else:
                f.write(
                    f"{gene_id} {gene.chr_id} {gene.start - 1} {gene.end - 1} {gene.strand}\n")
    return genefile



class LDAK_GBAT_JOB():
    def __init__(self, name, work_dir, input_vcf_file, input_gwas_df_h5, input_gff3_file):
        self.name = name
        self.work_dir = work_dir
        self.input_vcf_file = input_vcf_file
        self.input_gwas_df_h5 = input_gwas_df_h5
        self.input_gff3_file = input_gff3_file

    def prepare_data(self):
        """
        Prepare the data for the LDAK analysis
        1. Generate genetype file: self.genetype_file
        2. Generate summary statistics file: self.sum_stat_file
        3. Generate genefile: self.genefile
        """

        mkdir(self.work_dir)

        # get genetype file
        output_genetype_stem = self.input_vcf_file.replace(
            ".vcf.gz", ".ldak")
        if not have_file(output_genetype_stem + ".bed"):
            print("Generating genetype file")
            self.genetype_file = get_gentype_file(
                self.input_vcf_file, output_genetype_stem)
        else:
            print("Genetype file already exists")
            self.genetype_file = output_genetype_stem

        # get summary statistics file
        sum_stat_txt_file = self.input_gwas_df_h5.replace(
            ".h5", ".ldak.summaries")
        extract_file = self.input_gwas_df_h5.replace(".h5", ".ldak.extract")
        if not have_file(sum_stat_txt_file) or not have_file(extract_file):
            print("Generating summary statistics file")
            gwas_df = pd.read_hdf(self.input_gwas_df_h5)
            sample_number = len(get_sample_list_from_vcf(self.input_vcf_file))
            self.sum_stat_file, self.extract_file = get_summary_statistics_file(
                gwas_df, sample_number, sum_stat_txt_file, extract_file, chr_id_map_file=self.genetype_file + ".chr_id_map")
        else:
            print("Summary statistics file already exists")
            self.sum_stat_file = sum_stat_txt_file
            self.extract_file = extract_file

        # get genefile
        genefile = self.input_vcf_file.replace(".vcf.gz", ".ldak.genefile")
        if not have_file(genefile):
            print("Generating genefile")
            self.genefile = get_genefile(
                self.input_gff3_file, genefile, chr_id_map_file=self.genetype_file + ".chr_id_map")
        else:
            print("Genefile already exists")
            self.genefile = genefile

    def run_cut_genes(self):
        gene_buffer = 2000
        self.ldak_gbat = self.input_vcf_file.replace(".vcf.gz", ".ldak.gbat")

        if have_file(self.ldak_gbat):
            print("GBAT dir already exists")
        else:
            bash_script = f"{self.work_dir}/cut_genes.sh"
            with open(bash_script, 'w') as f:
                f.write(f"#!/bin/bash\n")
                f.write(f"source activate {LDAK_ENV}\n")
                f.write(f"{LDAK_PATH} --cut-genes {self.ldak_gbat} --by-chr YES --gene-buffer {gene_buffer} --bfile {self.genetype_file} --genefile {self.genefile} --allow-multi YES\n")
                f.write(f"source deactivate\n")
            
            cmd_run(f"bash {bash_script}", cwd=self.work_dir)

    def run(self):
        # prepare
        self.prepare_data()
        # run ldak
        ## cut genes
        self.run_cut_genes()


if __name__ == "__main__":

    from yxquantgene.analyses.gwas2 import GWAS_Job
    from yxquantgene.utils.vcf import concat_vcfs, quality_control_vcf, get_chr_id_list_from_vcf
    from yxutil import mkdir, cmd_run

    snp_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces_snp.vcf.gz"
    indel_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces_indel.vcf.gz"

    # merge snp and indel vcf
    concat_vcf_prefix = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces.uniq"
    concat_vcf_file = concat_vcfs([snp_vcf_file, indel_vcf_file], concat_vcf_prefix)

    # variant quality control
    min_maf = 0.01
    max_het = 0.1
    max_miss = 0.5

    vcf_qc_file = concat_vcf_file.replace(".vcf.gz", f".qc.miss{max_miss}.maf{min_maf}.het{max_het}.vcf.gz")
    vcf_qc_file = quality_control_vcf(concat_vcf_file, vcf_qc_file, min_maf=min_maf, max_het=max_het, max_miss=max_miss)

    # chrom rename
    chr_id_list = get_chr_id_list_from_vcf(vcf_qc_file)
    chr_id_map_file = vcf_qc_file.replace(".vcf.gz", ".chr_id_map")
    with open(chr_id_map_file, 'w') as f:
        for i, chr_id in enumerate(chr_id_list):
            f.write(f"{chr_id}\t{i+1}\n")
    qc_renamed_vcf_file = vcf_qc_file.replace(".vcf.gz", ".chr_rename.vcf.gz")
    cmd_run(f"bcftools annotate --rename-chrs {chr_id_map_file} {vcf_qc_file} -O z -o {qc_renamed_vcf_file}")
    cmd_run(f"tabix -f -p vcf {qc_renamed_vcf_file}")
    # bcftools annotate --rename-chrs /lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.chr_id_map.txt landraces.uniq.qc.miss0.5.maf0.01.het0.1.vcf.gz -O z -o landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_id.vcf.gz

    # plink
    bfile_prefix = qc_renamed_vcf_file.replace(".vcf.gz", "")
    cmd_run(f"plink --vcf {qc_renamed_vcf_file} --make-bed --out {bfile_prefix}")
    
    # magma annotation
    bfile_bim_file = bfile_prefix + ".bim"
    gene_loc_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.genefile"
    cmd_run(f"magma --annotate nonhuman window=2,1 --snp-loc {bfile_bim_file} --gene-loc {gene_loc_file} --out {bfile_prefix}")
    magma_annot_file = bfile_prefix + ".genes.annot"
    
    # write phenotype file
    bfile_fam_file = bfile_prefix + ".fam"
    pheno_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/pheno.fam"
    landraces_reseq_meta_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces_reseq_meta.csv"
    
    landraces_reseq_df = pd.read_csv(landraces_reseq_meta_file)
    bfile_fam_df = pd.read_csv(bfile_fam_file, sep='\s+', header=None)
    bfile_sample_list = bfile_fam_df[1].to_list()
    filtered_df = landraces_reseq_df[landraces_reseq_df['LIB'].isin(bfile_sample_list)]
    landraces_reseq_df = filtered_df.set_index('LIB').loc[bfile_sample_list].reset_index()

    pheno_name = 'ai'

    bfile_fam_df[2] = landraces_reseq_df[pheno_name]
    bfile_fam_df.iloc[:, :3].to_csv(pheno_file, sep='\t', header=False, index=False, na_rep='NA')
    
    # raw data gene 
    cmd_run(f"magma --bfile {bfile_prefix} --gene-annot {magma_annot_file} --pheno file={pheno_file} --out {bfile_prefix}.magme.results --burden 0.01")
    "magma --bfile landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename --pheno file=/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/pheno.fam --gene-annot landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename.genes.annot --out results_pheno4 --burden 0.01"

    # snp p value
    gemma_work_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename.gemma"
    gemma_phenotype_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/pheno.bimbam'
    kinship_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/population_structure/landraces_snp_pruned.win15000.miss0.5.maf0.01.het0.01.rq0.5.ibs.matrix"
    ref_genome_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    gemma_job = GWAS_Job(name='ai', work_dir=gemma_work_dir + "/ai", phenotype_file=gemma_phenotype_file,
                         genotype_file=qc_renamed_vcf_file, kinship_file=kinship_file, gwas_method='gemma', genome_fasta_file=ref_genome_file)
    gemma_job.run()
    gemma_job.gwas_df_h5

    gwas_df = pd.read_hdf(gemma_job.gwas_df_h5)
    p_value_file = gemma_job.gwas_df_h5.replace(".h5", ".p_value.txt")
    gwas_df[['ID', 'P_VALUE']].to_csv(p_value_file, sep='\t', index=False, header=False)
    
    "magma --bfile landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename --gene-annot landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename.genes.annot --pval landraces.uniq.qc.miss0.5.maf0.01.het0.1.chr_rename.gemma/ai/ai.gemma.gwas_df.p_value.txt N=443 --out results_ai_multi_adappermp --gene-model multi=snp-wise --gene-settings adap-permp"




    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces.uniq.vcf.gz"
    ref_genome_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    phenotype_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/pheno.bimbam'
    kinship_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/population_structure/landraces_snp_pruned.win15000.miss0.5.maf0.01.het0.01.rq0.5.ibs.matrix"

    work_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/ldak_test"
    mkdir(work_dir)

    gemma_job = GWAS_Job(name='ai', work_dir=work_dir + "/ai", phenotype_file=phenotype_file,
                         genotype_file=input_vcf_file, kinship_file=kinship_file, gwas_method='gemma', genome_fasta_file=ref_genome_file)
    gemma_job.run()
    gemma_job.gwas_df_h5

    # run ldak
    gff3_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3"

    ldak_job = LDAK_GBAT_JOB(name='ai', work_dir=work_dir + "/ai",
                             input_vcf_file=input_vcf_file, input_gwas_df_h5=gemma_job.gwas_df_h5, input_gff3_file=gff3_file)
    ldak_job.run()

    # run gbat
