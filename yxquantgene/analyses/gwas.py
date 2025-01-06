from yxutil import have_file, read_list_file
from yxseq import read_fasta_by_faidx
from statsmodels.sandbox.stats.multicomp import multipletests
from yxquantgene.plot.qqplot import quantile_quantile_plot
from yxquantgene.plot.manhattan import genomewide_gwas_plot as manhattan_plot
from yxquantgene.plot.geogenoplot import allele_corr_plot
import matplotlib.pyplot as plt
import pandas as pd
import pysam


class Variant():
    def __init__(self, chr_id=None, pos=None, ref=None, alt=None):
        self.chr_id = chr_id
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def load_from_gemma(self, gemma_result_row, vcf_file_path):
        self.chr_id = gemma_result_row['chr']
        self.pos = gemma_result_row['pos']
        self.ref = gemma_result_row['ref']
        self.alt = gemma_result_row['alt']
        self.minor = gemma_result_row['minor']
        self.major = gemma_result_row['major']
        self.n_miss = gemma_result_row['n_miss']
        self.maf = gemma_result_row['maf']
        self.beta = gemma_result_row['beta']
        self.alt_beta = self.beta if self.ref == self.major else -self.beta
        self.se = gemma_result_row['se']
        self.logl_H1 = gemma_result_row['logl_H1']
        self.l_remle = gemma_result_row['l_remle']
        self.l_mle = gemma_result_row['l_mle']
        self.pval = gemma_result_row['pval']
        self.fdr_bh = gemma_result_row['fdr_bh']
        self.get_vcf_rec(vcf_file_path)
        self.phase_variant_INFO()
        self.get_most_impact_snpeff_info()

        self.high_impact = self.most_impact_snpeff_info['annotation_impact'] == 'HIGH'
        self.annotation_impact = self.most_impact_snpeff_info['annotation_impact']
        self.in_gene = self.most_impact_snpeff_info['feature_type'] != 'intergenic_region'
        self.annotation = self.most_impact_snpeff_info['annotation']
        self.gene_id = self.most_impact_snpeff_info['gene_id']

    def get_donors(self, top_num=20, genotype='alt', reverse=True):
        """
        genotype: 'ref' or 'alt' or 'hetero'
        """
        if not hasattr(self, 'genotype_vs_phenotype_df'):
            raise ValueError("Please run get_genotype_vs_phenotype first")

        if genotype == 'ref':
            donor_df = self.genotype_vs_phenotype_df[self.genotype_vs_phenotype_df["genotype"] == "0|0"].sort_values(
                'value', ascending=(not reverse)).head(top_num)
        elif genotype == 'alt':
            donor_df = self.genotype_vs_phenotype_df[self.genotype_vs_phenotype_df["genotype"] == "1|1"].sort_values(
                'value', ascending=(not reverse)).head(top_num)
        elif genotype == 'hetero':
            donor_df = self.genotype_vs_phenotype_df[(self.genotype_vs_phenotype_df["genotype"] == "0|1") | (self.genotype_vs_phenotype_df["genotype"] == "1|0")].sort_values(
                'value', ascending=(not reverse)).head(top_num)

        return donor_df

    def plot(self, save_file=None, y_label='Value'):
        if not hasattr(self, 'genotype_vs_phenotype_df'):
            raise ValueError("Please run get_genotype_vs_phenotype first")

        ax, p0, p1 = allele_corr_plot(self.genotype_vs_phenotype_df, ref_allele=self.ref, alt_allele=self.alt,
                                      phased=False, trendline=True, y_label=y_label, title='Correlation', save_file=save_file)

        self.pearson_correlation = p0
        self.pearson_pvalue = p1

    def get_genotype_vs_phenotype(self, phenotype_file, acc_id_file):

        phenotype_list = read_list_file(phenotype_file)
        acc_id_list = read_list_file(acc_id_file)

        genotype_list = self.vcf_rec.samples

        genotype_vs_phenotype_df = pd.DataFrame(columns=[
            'acc_id', 'genotype', 'value'])

        for i in range(len(acc_id_list)):
            if phenotype_list[i] == "NA":
                continue
            genotype_vs_phenotype_df.loc[i] = [acc_id_list[i], "|".join(
                ["0" if j == self.ref else "1" for j in genotype_list[acc_id_list[i]].alleles]), float(phenotype_list[i])]

        self.genotype_vs_phenotype_df = genotype_vs_phenotype_df

    def get_vcf_rec(self, vcf_file_path):
        chrom = self.chr_id
        pos = self.pos
        ref = self.ref
        alt = self.alt

        with pysam.VariantFile(vcf_file_path) as vcf:
            for rec in vcf.fetch(chrom, pos-1, pos):
                if rec.pos == pos:
                    if ref is not None and alt is not None:
                        if rec.ref == ref and rec.alts[0] == alt:
                            self.vcf_rec = rec
                            return

    def phase_variant_INFO(self):
        """
        ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
        """
        info_dict = {}
        num = 0
        for info_str in self.vcf_rec.info['ANN']:
            data_list = info_str.split("|")
            info_dict[num] = {
                'allele': data_list[0],
                'annotation': data_list[1].split("&"),
                # HIGH, MODERATE, LOW, MODIFIER, choose HIGH for most cases
                'annotation_impact': data_list[2],
                'gene_name': data_list[3],
                'gene_id': data_list[4],
                'feature_type': data_list[5],
                'feature_id': data_list[6],
                'transcript_biotype': data_list[7],
                'rank': data_list[8],
                'hgvs_c': data_list[9],
                'hgvs_p': data_list[10],
                'cdna_pos_cdna_length': data_list[11],
                'cds_pos_cds_length': data_list[12],
                'aa_pos_aa_length': data_list[13],
                'distance': data_list[14],
                'errors_warnings_info': data_list[15]
            }
            num += 1
        self.snpeff_info = info_dict

    def get_most_impact_snpeff_info(self, impact_sort_list=['HIGH', 'MODERATE', 'LOW', 'MODIFIER']):
        self.most_impact_snpeff_info = None
        for impact in impact_sort_list:
            for info in self.snpeff_info.values():
                if info['annotation_impact'] == impact:
                    self.most_impact_snpeff_info = info
                    return


class Gemma_Job():
    def __init__(self, name=None, work_dir=None, phenotype_file=None, genotype_file=None, kinship_file=None, raw_vcf_file=None, acc_id_file=None, genome_fasta_file=None, genome_gff_file=None):
        self.name = name
        self.phenotype_file = phenotype_file
        self.genotype_file = genotype_file
        self.kinship_file = kinship_file
        self.work_dir = work_dir
        self.raw_vcf_file = raw_vcf_file
        self.acc_id_file = acc_id_file
        self.genome_fasta_file = genome_fasta_file
        self.genome_gff_file = genome_gff_file

        if self.acc_id_file:
            self.acc_id_list = read_list_file(self.acc_id_file)
        else:
            self.acc_id_list = None

    def print_cmd_file(self, cmd_file):
        cms_string = "gemma -g %s -k %s -lmm 4 -miss 0.1 -p %s -o %s\n" % (
            self.genotype_file, self.kinship_file, self.phenotype_file, self.name)
        with open(cmd_file, "a") as f:
            f.write("source activate gemma\n")
            f.write("cd %s\n" % self.work_dir)
            f.write(cms_string)

    def save(self, save_file=None):
        if save_file:
            self.assoc_df.to_csv(save_file, sep=",", index=False)
        else:
            save_file = self.work_dir + "/%s.assoc.csv" % self.name
            self.assoc_df.to_csv(save_file, sep=",", index=False)

    def load_from_csv(self, csv_file=None):
        if csv_file is None:
            csv_file = self.work_dir + "/%s.assoc.csv" % self.name
        self.assoc_df = pd.read_csv(csv_file)

    def load_from_gemma(self, assoc_file=None, log_file=None):
        if assoc_file:
            self.output_assoc_file = assoc_file
        if log_file:
            self.gemma_log_file = log_file

        if assoc_file is None:
            output_dir = self.work_dir + "/output"
            self.output_assoc_file = output_dir + "/%s.assoc.txt" % self.name
            self.gemma_log_file = output_dir + "/%s.log.txt" % self.name

        if not have_file(self.output_assoc_file):
            print("Error: %s not found" % self.output_assoc_file)
        if not have_file(self.gemma_log_file):
            print("Error: %s not found" % self.gemma_log_file)

        self.gemma_assoc = pd.read_csv(self.output_assoc_file, sep="\t")

        self.assoc_df = pd.DataFrame()

        self.assoc_df['chr'] = self.gemma_assoc['rs'].map(
            lambda x: x.split("_")[0])
        self.assoc_df['pos'] = self.gemma_assoc['rs'].map(
            lambda x: int(x.split("_")[1]))
        self.assoc_df['ref'] = self.gemma_assoc['rs'].map(
            lambda x: x.split("_")[2].split("/")[0])
        self.assoc_df['alt'] = self.gemma_assoc['rs'].map(
            lambda x: x.split("_")[2].split("/")[1])
        self.assoc_df['minor'] = self.gemma_assoc['allele1']
        self.assoc_df['major'] = self.gemma_assoc['allele0']
        self.assoc_df['n_miss'] = self.gemma_assoc['n_miss']
        self.assoc_df['maf'] = self.gemma_assoc['af']
        self.assoc_df['beta'] = self.gemma_assoc['beta']
        self.assoc_df['se'] = self.gemma_assoc['se']
        self.assoc_df['logl_H1'] = self.gemma_assoc['logl_H1']
        self.assoc_df['l_remle'] = self.gemma_assoc['l_remle']
        self.assoc_df['l_mle'] = self.gemma_assoc['l_mle']
        # likelihood ratio test
        self.assoc_df['pval'] = self.gemma_assoc['p_lrt']

    def multiple_test(self, method='fdr_bh'):
        self.assoc_df[method] = multipletests(
            self.assoc_df['pval'], method=method)[1]

    def browse_variant(self, sort_by='pval'):
        self.multiple_test()
        assoc_df = self.assoc_df.sort_values(sort_by)

        for i, row in assoc_df.iterrows():
            var = Variant()
            var.load_from_gemma(row, self.raw_vcf_file)
            yield var

    def qqplot(self, plot_keys=['pval', 'bonferroni', 'fdr_bh'], subtitles=['non-adjusted P-value', 'Bonferroni', 'FDR'], save_file=None):
        for k in plot_keys:
            if k not in self.assoc_df.columns:
                print(
                    "Error: %s not in assoc_df dataframe, please run multiple_test first" % k)
                return

        fig, axs = plt.subplots(1, len(plot_keys), figsize=(12, 6))

        for i, k in enumerate(plot_keys):
            quantile_quantile_plot(self.assoc_df[k], subtitles[i], ax=axs[i])

        plt.tight_layout()

        if save_file:
            fig.savefig(save_file, format='png', facecolor='none', dpi=300,
                        edgecolor='none', bbox_inches='tight')

        plt.show()

    def manhattan_plot(self, multitest_method='fdr_bh', threshold_qval=0.05, min_plot_chr_length=1e6, save_file=None):
        if not self.genome_fasta_file:
            print("Error: genome_fasta_file not set")
            return

        genome_dict = read_fasta_by_faidx(self.genome_fasta_file)
        chr_length_dict = {k: v.len() for k, v in genome_dict.items()}
        chr_list = list(
            [i for i in chr_length_dict if chr_length_dict[i] > min_plot_chr_length])

        # 取出job.assoc_df 中的['chr', 'pos', 'pval', 'fdr_bh']列，赋值给manhattan_df
        manhattan_df = self.assoc_df[['chr', 'pos', 'pval', multitest_method]]
        manhattan_df.columns = ['chr', 'pos', 'pval', 'qval']

        fig, ax = plt.subplots(figsize=(16, 9))
        manhattan_plot(manhattan_df, chr_list, chr_length_dict,
                       threshold_fdr=threshold_qval, ax=ax)
        if save_file:
            fig.savefig(save_file, format='png', facecolor='none', dpi=300,
                        edgecolor='none', bbox_inches='tight')
        plt.show()


if __name__ == "__main__":
    job_name = "jesse_ai_snp"
    phenotype_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/3.envGWAS/jesse/ai/pheno.bimbam"
    genotype_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/3.envGWAS/jesse/snp.bimbam"
    kinship_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/3.envGWAS/jesse/kinship.ibs.txt"
    gemma_work_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/3.envGWAS/jesse/ai/snp"
    raw_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/Sorghum_d8.noduplicates.allChr.snp._markernamesadded_imputed_snpeff.vcf.gz"
    acc_id_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/3.envGWAS/jesse/acc_id.txt"
    genome_fasta_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    genome_gff_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3"

    # load job
    job = Gemma_Job(job_name, gemma_work_dir)
    job.load_from_gemma()
    print(job.assoc_df.head())

    # multiple test
    job.multiple_test('bonferroni')
    job.multiple_test('fdr_bh')
    print(job.assoc_df.head())
    qqplot_file = job.work_dir + "/qqplot.png"
    job.qqplot(save_file=qqplot_file)

    # save job
    job.save()

    # manhattan plot
    job.genome_fasta_file = genome_fasta_file
    manhattan_file = job.work_dir + "/manhattan.png"
    job.manhattan_plot(save_file=manhattan_file)

    # get variant info
    job = Gemma_Job(job_name, gemma_work_dir, raw_vcf_file=raw_vcf_file,
                    acc_id_file=acc_id_file, phenotype_file=phenotype_file)
    job.load_from_gemma()

    chr_id = "Chr08"
    pos = 16762784
    ref = "T"
    alt = "C"

    gemma_result_row = job.assoc_df.loc[(job.assoc_df['chr'] == chr_id) & (
        job.assoc_df['pos'] == pos) & (job.assoc_df['ref'] == ref) & (job.assoc_df['alt'] == alt)].iloc[0]

    var = Variant()
    var.load_from_gemma(gemma_result_row, job.raw_vcf_file)
    var.get_genotype_vs_phenotype(phenotype_file, acc_id_file)
    var.plot()
    var.get_candidate_donor(positive=True, top_num=20)

    print(var.ref, var.alt, var.minor, var.major)
    print(var.beta)
    print(var.pearson_correlation)
