from yxutil import cmd_run, mkdir, have_file
from yxquantgene.utils.vcf import get_chr_id_list_from_vcf, get_sample_list_from_vcf
from yxquantgene.utils.utils import get_gene_range_df
from statsmodels.sandbox.stats.multicomp import multipletests
from yxseq import read_gff_file
import pandas as pd


def results_parser(magma_output_genes_out, input_gff3_file, model='snp_pvalue'):
    gene_range_df = get_gene_range_df(input_gff3_file)
    gene_range_df = gene_range_df.rename(columns={'GENE': 'ID'})

    if model == 'snp_pvalue':
        magma_output_genes_df = pd.read_csv(
            magma_output_genes_out, skiprows=1, delim_whitespace=True)
        magma_output_genes_df = magma_output_genes_df.rename(
            columns={'GENE': 'ID', 'PERMP_MULTI': 'P_VALUE'})
        gwas_df = magma_output_genes_df[['ID', 'NSNPS', 'P_VALUE']]
        gwas_df['FDR'] = multipletests(
            gwas_df['P_VALUE'], method='fdr_bh')[1]
        gwas_df = pd.merge(gwas_df, gene_range_df, on='ID', how='left')

    return gwas_df


def get_genefile(input_gff3_file, genefile, chr_id_map_file=None):
    if chr_id_map_file:
        chr_id_map_dict = {}
        with open(chr_id_map_file) as f:
            for line in f:
                old_name, new_name = line.strip().split()
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


class MAGMA_JOB():
    def __init__(self, name, work_dir, input_vcf_file, gff_file, site_gwas_df_h5=None, phenotype_file=None):
        self.name = name
        self.work_dir = work_dir
        self.input_vcf_file = input_vcf_file
        self.gff_file = gff_file
        self.phenotype_file = phenotype_file
        self.site_gwas_df_h5 = site_gwas_df_h5

    def prepare(self):
        chr_id_map_file = self.input_vcf_file.replace(".vcf.gz", ".chr_id_map")
        renamed_vcf_file = self.input_vcf_file.replace(
            ".vcf.gz", ".chr_rename.vcf.gz")
        bfile_prefix = renamed_vcf_file.replace(".vcf.gz", "")

        if not (have_file(bfile_prefix + ".bed") and have_file(chr_id_map_file) and have_file(renamed_vcf_file)):
            # chrom rename
            chr_id_list = get_chr_id_list_from_vcf(self.input_vcf_file)

            with open(chr_id_map_file, 'w') as f:
                for i, chr_id in enumerate(chr_id_list):
                    f.write(f"{chr_id}\t{i+1}\n")

            cmd_run(
                f"bcftools annotate --rename-chrs {chr_id_map_file} {self.input_vcf_file} -O z -o {renamed_vcf_file}")
            cmd_run(f"tabix -f -p vcf {renamed_vcf_file}")

            # plink
            cmd_run(
                f"plink --vcf {renamed_vcf_file} --make-bed --out {bfile_prefix}")

        self.chr_id_map_file = chr_id_map_file
        self.renamed_vcf_file = renamed_vcf_file
        self.bfile_prefix = bfile_prefix
        self.sample_num = len(get_sample_list_from_vcf(self.renamed_vcf_file))

    def make_p_value_file(self):
        site_gwas_df = pd.read_hdf(self.site_gwas_df_h5)
        p_value_file = self.site_gwas_df_h5.replace(".h5", ".p_value.txt")
        site_gwas_df[['ID', 'P_VALUE']].to_csv(
            p_value_file, sep='\t', index=False, header=False)
        self.p_value_file = p_value_file

    def annotate(self):
        bfile_bim_file = self.bfile_prefix + ".bim"
        self.gene_annot_file = self.bfile_prefix + ".genes.annot"

        if not have_file(self.gene_annot_file):
            gene_loc_file = self.bfile_prefix + ".ref_gene_loc"
            get_genefile(self.gff_file, gene_loc_file, self.chr_id_map_file)
            cmd_run(
                f"magma --annotate nonhuman window=2,1 --snp-loc {bfile_bim_file} --gene-loc {gene_loc_file} --out {self.bfile_prefix}")

    def gene_analysis(self, model='snp_pvalue'):
        if model == 'snp_pvalue':
            self.make_p_value_file()
            self.magma_output_prefix = f"{self.output_dir}/magma.snp_pvalue.results"
            if not have_file(self.magma_output_prefix + ".genes.out"):
                cmd_run(f"magma --bfile {self.bfile_prefix} --gene-annot {self.gene_annot_file} --pval {self.p_value_file} N={self.sample_num} --out {self.magma_output_prefix} --gene-model multi=snp-wise --gene-settings adap-permp")

    def run(self, model='snp_pvalue'):
        self.output_dir = f"{self.work_dir}/{self.name}.magma"
        mkdir(self.output_dir)

        self.prepare()
        self.annotate()
        self.gene_analysis(model=model)

        # 解析结果
        gwas_df_h5 = f'{self.work_dir}/{self.name}.magma.gwas_df.h5'
        if not have_file(gwas_df_h5):
            gwas_df = results_parser(
                f'{self.magma_output_prefix}.genes.out', self.gff_file, model=model)
            gwas_df.to_hdf(
                gwas_df_h5, key='df', mode='w')

        self.gwas_df_h5 = gwas_df_h5

    def __str__(self):
        return f"MAGMA JOB: {self.name}"


if __name__ == "__main__":
    # input vcf file
    # A VCF file with genotypes, sample order should be the same as other files
    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/landraces_indel.vcf.gz"
    # gff_file
    # A standard GFF3 file
    gff_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3"
    # phenotype_file
    # A single column CSV file without header, values as phenotype, sample order should be the same as vcf_file
    # phenotype_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/pheno.bimbam"
    # site_gwas_df_h5
    # A HDF5 file with key 'df' and columns:
    # For site gwas, gwas_df should have columns: ID, CHROM, POS, REF, ALT, BETA, P_VALUE, FDR
    site_gwas_df_h5 = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/ai.gemma.gwas_df.h5"

    # running
    # MAGMA_JOB
    name = "ai"
    work_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test"
    mkdir(work_dir)

    magma_job = MAGMA_JOB(name=name, work_dir=work_dir, input_vcf_file=input_vcf_file,
                          gff_file=gff_file, site_gwas_df_h5=site_gwas_df_h5)
    magma_job.run()
