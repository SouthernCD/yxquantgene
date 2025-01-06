from statsmodels.sandbox.stats.multicomp import multipletests
from time import sleep
from yxquantgene.plot.manhattan import genomewide_gwas_plot as manhattan_plot
from yxquantgene.utils.vcf import convert_vcf_to_bimbam, build_var_stat_table, read_hdf_with_retry, get_chr_id_list_from_vcf
from yxquantgene.utils.utils import get_chr_len_df, get_gene_range_df
from yxquantgene.metrics.wza import gene_wza_pipeline
from yxquantgene.analyses.gwas_wrapper.gemma_wrapper import GEMMA_JOB
from yxquantgene.analyses.gwas_wrapper.magma_wrapper import MAGMA_JOB
from yxutil import have_file, read_list_file, mkdir, cmd_run, multiprocess_running, pickle_load, pickle_dump
import matplotlib.pyplot as plt
import pandas as pd
import random
import numpy as np

# normal GWAS


class GWAS_Job():
    def __init__(self, name=None, work_dir=None, **kwargs):
        self.name = name
        self.phenotype_file = kwargs.get('phenotype_file', None)
        self.genotype_file = kwargs.get('genotype_file', None)
        self.kinship_file = kwargs.get('kinship_file', None)
        self.gwas_method = kwargs.get('gwas_method', 'gemma')
        self.work_dir = work_dir
        self.acc_id_file = kwargs.get('acc_id_file', None)
        self.genome_fasta_file = kwargs.get('genome_fasta_file', None)
        self.genome_gff_file = kwargs.get('genome_gff_file', None)

        # 创建工作目录
        mkdir(self.work_dir)

        # 读取样品ID
        if self.acc_id_file:
            self.acc_id_list = read_list_file(self.acc_id_file)
        else:
            self.acc_id_list = None

        # # 计算基因组染色体长度
        # if self.genome_fasta_file:
        #     self.chr_len_df = get_chr_len_df(self.genome_fasta_file)
        # else:
        #     self.chr_len_df = None

        # # 计算基因范围
        # if self.genome_gff_file:
        #     self.gene_range_df = get_gene_range_df(self.genome_gff_file)
        # else:
        #     self.gene_range_df = None

    def run(self, **kwargs):
        if self.gwas_method == 'gemma':
            retry = kwargs.get('retry', 3)
            gemma_job = GEMMA_JOB(name=self.name, work_dir=self.work_dir, phenotype_file=self.phenotype_file,
                                  input_vcf_file=self.genotype_file, kinship_file=self.kinship_file)
            gemma_job.run(retry=retry)
            self.gwas_df_h5 = gemma_job.gwas_df_h5
            return self.gwas_df_h5
        else:
            raise ValueError(f'Unsupported GWAS method: {self.gwas_method}')

    def manhattan_plot(self, genome_fasta_file):
        if not hasattr(self, 'gwas_df_h5'):
            raise ValueError(
                'Please run GWAS first before plotting Manhattan plot.')

        # 计算基因组染色体长度
        chr_len_df = get_chr_len_df(genome_fasta_file)

        # 绘制Manhattan图
        self.gwas_manhattan_plot_file = self.gwas_df_h5.replace(
            '.h5', '.manhattan.png')

        if not have_file(self.gwas_manhattan_plot_file):
            gwas_df = pd.read_hdf(self.gwas_df_h5, key='df')
            manhattan_df = gwas_df[['CHROM', 'POS', 'P_VALUE', 'FDR']]

            fig, ax = plt.subplots(figsize=(16, 9))
            manhattan_plot(manhattan_df, chr_len_df, threshold_fdr=0.05, ax=ax)
            ax.set_title(f'{self.name} GWAS Manhattan Plot')
            fig.savefig(self.gwas_manhattan_plot_file, format='png', facecolor='none', dpi=300,
                        edgecolor='none', bbox_inches='tight')

        return self.gwas_manhattan_plot_file


# site to gene level GWAS


class GeneGWAS_JOB():
    def __init__(self, name=None, work_dir=None, **kwargs):
        self.name = name
        self.work_dir = work_dir
        self.genome_gff_file = kwargs.get('genome_gff_file', None)
        self.gene_flank = kwargs.get('gene_flank', 2000)
        self.site2gene_method = kwargs.get('site2gene_method', 'wza')
        self.site_gwas_df_h5 = kwargs.get('site_gwas_df_h5', None)
        self.phenotype_file = kwargs.get('phenotype_file', None)
        self.genotype_file = kwargs.get('genotype_file', None)
        self.kinship_file = kwargs.get('kinship_file', None)
        self.acc_id_file = kwargs.get('acc_id_file', None)
        self.gwas_method = kwargs.get('gwas_method', 'gemma')
        self.n_jobs = kwargs.get('n_jobs', 1)

        # 创建工作目录
        mkdir(self.work_dir)

    def run(self):
        self.gene_gwas_df_h5 = f'{self.work_dir}/{self.name}.{self.site2gene_method}.gwas_df.h5'
        if have_file(self.gene_gwas_df_h5):
            print(f'Gene level GWAS already exists.')
            return self.gene_gwas_df_h5
        else:
            # load site level GWAS results
            if self.site_gwas_df_h5 is None or not have_file(self.site_gwas_df_h5):
                gwas_job = GWAS_Job(name=self.name, work_dir=self.work_dir, phenotype_file=self.phenotype_file,
                                    genotype_file=self.genotype_file, kinship_file=self.kinship_file, acc_id_file=self.acc_id_file)
                gwas_job.run(gwas_method=self.gwas_method)

                self.site_gwas_df_h5 = gwas_job.gwas_df_h5
                
            if self.site2gene_method == 'wza':
                self.site_gwas_df = pd.read_hdf(self.site_gwas_df_h5, key='df')
                self.gene_gwas_df = self.site2gene_wza()
            elif self.site2gene_method == 'magma':
                magma_job = MAGMA_JOB(name=self.name, work_dir=self.work_dir, input_vcf_file=self.genotype_file,
                                    gff_file=self.genome_gff_file, site_gwas_df_h5=self.site_gwas_df_h5)
                magma_job.run()
                self.gene_gwas_df = pd.read_hdf(magma_job.gwas_df_h5, key='df')
            else:
                raise ValueError(
                    f'Unsupported site2gene method: {self.site2gene_method}')

            self.gene_gwas_df.to_hdf(self.gene_gwas_df_h5, key='df', mode='w')

            return self.gene_gwas_df_h5

    def site2gene_wza(self):
        """
        site_gwas_df: pandas.DataFrame
            A DataFrame with columns: 'CHROM', 'POS', 'P_VALUE', 'MAF'
        gene_range_df: dict
            A DataFrame with columns: 'GENE', 'CHROM', 'START', 'END'
        """
        # load gene region
        self.gene_range_df = get_gene_range_df(self.genome_gff_file)

        gene_df = pd.DataFrame()
        gene_df['GENE'] = self.gene_range_df['GENE']
        gene_df['CHROM'] = self.gene_range_df['CHROM']
        gene_df['START'] = self.gene_range_df['START'] - self.gene_flank
        gene_df['END'] = self.gene_range_df['END'] + self.gene_flank

        gene_df = gene_wza_pipeline(
            self.site_gwas_df, gene_df, n_jobs=self.n_jobs)

        return gene_df

    def manhattan_plot(self, genome_fasta_file):
        if not hasattr(self, 'gene_gwas_df_h5'):
            raise ValueError(
                'Please run GWAS first before plotting Manhattan plot.')

        # 计算基因组染色体长度
        chr_len_df = get_chr_len_df(genome_fasta_file)

        # 绘制Manhattan图
        self.gwas_manhattan_plot_file = self.gene_gwas_df_h5.replace(
            '.h5', '.manhattan.png')

        if not have_file(self.gwas_manhattan_plot_file):
            gene_gwas_df = pd.read_hdf(self.gene_gwas_df_h5, key='df')
            manhattan_df = gene_gwas_df[['CHROM', 'P_VALUE', 'FDR']]
            manhattan_df['POS'] = (
                gene_gwas_df['START'] + gene_gwas_df['END']) // 2

            fig, ax = plt.subplots(figsize=(16, 9))
            manhattan_plot(manhattan_df, chr_len_df, threshold_fdr=0.05, ax=ax)
            ax.set_title(f'{self.name} GWAS Manhattan Plot')
            fig.savefig(self.gwas_manhattan_plot_file, format='png', facecolor='none', dpi=300,
                        edgecolor='none', bbox_inches='tight')

        return self.gwas_manhattan_plot_file


class PermGWAS_JOB():
    def __init__(self, name=None, work_dir=None, raw_job=None, **job_kwargs):
        self.name = name
        self.work_dir = work_dir
        self.job_kwargs = job_kwargs
        self.raw_job = raw_job

        self.phenotype_file = job_kwargs.get('phenotype_file', None)
        self.genotype_file = job_kwargs.get('genotype_file', None)
        self.kinship_file = job_kwargs.get('kinship_file', None)
        self.acc_id_file = job_kwargs.get('acc_id_file', None)
        self.genome_fasta_file = job_kwargs.get('genome_fasta_file', None)
        self.genome_gff_file = job_kwargs.get('genome_gff_file', None)

        # 创建工作目录
        mkdir(self.work_dir)

    def run(self, n_perm=1000, n_jobs=10):
        """
        进行GWAS permutation
        """
        print(f'Running GWAS permutation...')
        print(f'Number of permutation: {n_perm}')

        if n_jobs > 30:
            n_jobs = 30
            print(f'Number of jobs is limited to 30.')

        # 创建pheno_perm文件目录
        perm_root_dir = f'{self.work_dir}/pheno_perm'
        mkdir(perm_root_dir)

        # 创建pheno_perm文件
        pheno_value_list = read_list_file(self.phenotype_file)

        args_dict = {}
        for i in range(n_perm):
            perm_dir = f'{perm_root_dir}/perm_{i}'
            mkdir(perm_dir)
            pheno_perm_file = f'{perm_dir}/pheno_perm.bimbam'
            if not have_file(pheno_perm_file):
                random.shuffle(pheno_value_list)
                with open(pheno_perm_file, 'w') as f:
                    for pheno_value in pheno_value_list:
                        f.write(f'{pheno_value}\n')
            permute_job = self.raw_job(
                name=f'perm_{i}', work_dir=perm_dir, **self.job_kwargs)
            permute_job.phenotype_file = pheno_perm_file

            args_dict[f'perm_{i}'] = (permute_job,)

        mlt_dict = multiprocess_running(permute_job_run, args_dict, n_jobs)

        perm_job_results_dict = {i: mlt_dict[i]['output'] for i in mlt_dict}

        self.perm_job_results_dict = perm_job_results_dict

    def permutation_test(self, ori_gwas_df_h5):
        results_df_csv = ori_gwas_df_h5.replace('.h5', '.perm_results.csv')

        if have_file(results_df_csv):
            return results_df_csv
        else:
            results_df = pd.DataFrame(
                columns=['THRESHOLD', 'OBSERVED', 'PERMUTE', 'TYPE', 'P_VALUE'])

            # site level permutation test
            gwas_df = pd.read_hdf(ori_gwas_df_h5, key='df')

            for i in range(1, 51, 1):
                threshold = i / 100
                for var in ['P_VALUE', 'FDR']:
                    observed_sig_site_count = len(
                        gwas_df[gwas_df[var] < threshold])
                    results_df.loc[len(results_df)] = [
                        threshold, observed_sig_site_count, "", var, 0]

            for pj_id in self.perm_job_results_dict:
                print(f'Processing permutation {pj_id}...')
                permute_assoc_df = pd.read_hdf(
                    self.perm_job_results_dict[pj_id], key='df')

                for i, r in results_df.iterrows():
                    threshold = r['THRESHOLD']
                    var = r['TYPE']
                    if len(r['PERMUTE']) == 0:
                        permute_sig_site_counts = []
                    else:
                        permute_sig_site_counts = [
                            int(j) for j in r['PERMUTE'].split(',')]

                    sig_sites = permute_assoc_df[permute_assoc_df[var]
                                                 < threshold].index
                    permute_sig_site_counts.append(len(sig_sites))
                    results_df.loc[i, 'PERMUTE'] = ','.join(
                        [str(j) for j in permute_sig_site_counts])

            results_df['P_VALUE'] = results_df.apply(
                lambda x: calc_permute_p_value(x['OBSERVED'], [int(i) for i in x['PERMUTE'].split(',')]), axis=1)

            results_df.to_csv(results_df_csv, index=False)

            return results_df_csv


def permute_job_run(job):
    return job.run()


def calc_permute_p_value(observed, permute):
    return sum([1 for i in permute if i >= observed]) / len(permute)


if __name__ == "__main__":
    # Running GWAS with GEMMA
    from yxquantgene import GWAS_Job
    work_dir = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test'

    genotype_file = f'{work_dir}/landraces_indel.vcf.gz'
    phenotype_file = f'{work_dir}/pheno.bimbam'
    kinship_file = f'{work_dir}/kinship.matrix.csv'

    name = 'ai'
    gwas_method = 'gemma'

    gwas_job = GWAS_Job(name=name, work_dir=work_dir, phenotype_file=phenotype_file,
                        genotype_file=genotype_file, kinship_file=kinship_file, gwas_method=gwas_method)
    gwas_df = gwas_job.run()
    gwas_job.gwas_df_h5

    genome_fasta_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    gwas_job.manhattan_plot(genome_fasta_file)

    gwas_job.gwas_manhattan_plot_file

    # from site to gene level GWAS
    from yxquantgene import GeneGWAS_JOB
    genome_gff_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3"
    # gene_flank = 2000
    site2gene_method = 'magma'

    gene_gwas_job = GeneGWAS_JOB(name=name, work_dir=work_dir, site2gene_method=site2gene_method, genome_gff_file=genome_gff_file,
                                 gwas_method=gwas_method, phenotype_file=phenotype_file, genotype_file=genotype_file, kinship_file=kinship_file)
    gene_gwas_df = gene_gwas_job.run()
    gene_gwas_job.gene_gwas_df_h5

    gene_gwas_job.manhattan_plot(genome_fasta_file)
    gene_gwas_job.gwas_manhattan_plot_file

    # permutation test
    from yxquantgene import PermGWAS_JOB

    name = 'ai_site_gwas_perm'
    work_dir = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test'
    genome_fasta_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    genome_gff_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3"
    
    genotype_file = f'{work_dir}/landraces_indel.vcf.gz'
    phenotype_file = f'{work_dir}/pheno.bimbam'
    kinship_file = f'{work_dir}/kinship.matrix.csv'

    gwas_method = 'gemma'

    permGWAS_job = PermGWAS_JOB(name=name, work_dir=work_dir, raw_job=GWAS_Job, gwas_method=gwas_method, phenotype_file=phenotype_file,
                                genotype_file=genotype_file, kinship_file=kinship_file, acc_id_file=None, genome_fasta_file=genome_fasta_file, genome_gff_file=genome_gff_file)
    permGWAS_job.run(n_perm=100, n_jobs=30)

    ori_gwas_df_h5 = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/ai.gemma.gwas_df.h5"
    perm_results_csv = permGWAS_job.permutation_test(ori_gwas_df_h5)

    # permutation test
    name = 'ai_gene_gwas_perm'
    gwas_method = 'gemma'
    site2gene_method = 'magma'

    permGWAS_job = PermGWAS_JOB(name=name, work_dir=work_dir, raw_job=GeneGWAS_JOB, site2gene_method=site2gene_method, gwas_method=gwas_method, phenotype_file=phenotype_file,
                                genotype_file=genotype_file, kinship_file=kinship_file, acc_id_file=None, genome_fasta_file=genome_fasta_file, genome_gff_file=genome_gff_file)
    permGWAS_job.run(n_perm=100, n_jobs=30)

    ori_gwas_df_h5 = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/ai.magma.gwas_df.h5"
    perm_results_csv = permGWAS_job.permutation_test(ori_gwas_df_h5)
