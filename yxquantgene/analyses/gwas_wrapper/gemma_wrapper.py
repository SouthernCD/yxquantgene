import pandas as pd
from yxutil import cmd_run, mkdir, have_file
from yxquantgene.utils.vcf import convert_vcf_to_bimbam, build_var_stat_table, read_hdf_with_retry, get_chr_id_list_from_vcf
from time import sleep
from statsmodels.sandbox.stats.multicomp import multipletests


def gemma_result_parser(gemma_output_assoc_file, var_stat_h5, chr_list):
    """
    从GEMMA的输出文件中解析出结果，并与var_stat_h5中的信息进行合并
    """

    print(f'Parsing {gemma_output_assoc_file}...')
    print(f'Loading gemma_output_assoc_file...')
    gemma_output_assoc_df = pd.read_csv(gemma_output_assoc_file, sep='\t')
    gemma_output_assoc_df = gemma_output_assoc_df[[
        'rs', 'beta', 'se', 'af', 'p_lrt']]
    gemma_output_assoc_df = gemma_output_assoc_df.rename(
        columns={'rs': 'ID', 'beta': 'BETA', 'se': 'SE', 'af': 'MAF', 'p_lrt': 'P_VALUE'})
    print(f'Calculating FDR...')
    gemma_output_assoc_df['FDR'] = multipletests(
        gemma_output_assoc_df['P_VALUE'], method='fdr_bh')[1]

    if var_stat_h5 is None:
        return gemma_output_assoc_df

    print(f'Loading var_stat_h5...')
    chr_pos_df_list = []
    for chr_name in chr_list:
        # print(f'Loading {chr_name}...')
        chr_pos_df = read_hdf_with_retry(var_stat_h5, chr_name)
        if len(chr_pos_df) == 0:
            continue
        chr_pos_df = chr_pos_df[['ID', 'CHROM', 'POS', 'MAF', 'REF', 'ALT']]
        chr_pos_df_list.append(chr_pos_df)

    print(f'Merging...')
    chr_pos_df = pd.concat(chr_pos_df_list)
    gemma_output_assoc_df = gemma_output_assoc_df.drop(columns=['MAF'])
    merged_df = pd.merge(gemma_output_assoc_df,
                         chr_pos_df, on='ID', how='left')

    return merged_df


class GEMMA_JOB():
    def __init__(self, name, work_dir, input_vcf_file, phenotype_file, kinship_file):
        self.name = name
        self.work_dir = work_dir
        self.input_vcf_file = input_vcf_file
        self.kinship_file = kinship_file
        self.phenotype_file = phenotype_file

    def prepare(self):
        var_stat_h5 = f'{self.input_vcf_file}.var_stat.h5'
        if not have_file(var_stat_h5):
            print(f'Building variant stat table...')
            build_var_stat_table(self.input_vcf_file, var_stat_h5)
        self.bimbam_file = f'{self.input_vcf_file}'.replace(
            '.vcf.gz', '.bimbam')
        if not have_file(self.bimbam_file):
            print(f'Converting VCF to BIMBAM...')
            convert_vcf_to_bimbam(self.input_vcf_file, self.bimbam_file)
        self.genome_chr_list = get_chr_id_list_from_vcf(self.input_vcf_file)

    def run(self, retry=3):
        self.prepare()

        # 创建输出目录
        output_dir = f'{self.work_dir}/{self.name}.gemma'
        mkdir(output_dir)

        # 运行GEMMA
        gemma_output_assoc_file = f'{output_dir}/output/gemma_output.assoc.txt'
        if have_file(gemma_output_assoc_file):
            print(
                f'GEMMA GWAS result already exists.')
        else:
            gemma_cmd_file = f'{output_dir}/gemma_cmd.sh'
            with open(gemma_cmd_file, 'w') as f:
                f.write(f'cd {output_dir}\n')
                f.write("source activate gemma\n")
                f.write(
                    f"gemma -g {self.bimbam_file} -k {self.kinship_file} -lmm 4 -miss 0.1 -p {self.phenotype_file} -o gemma_output\n")

            cmd_run(f'bash {gemma_cmd_file}', cwd=output_dir)

            retry_count = 0
            if not have_file(gemma_output_assoc_file):
                retry_count += 1
                if retry_count < retry:
                    print(
                        f'GEMMA GWAS failed, retrying {retry_count}...')
                    sleep(30)
                    cmd_run(f'bash {gemma_cmd_file}', cwd=output_dir)
                else:
                    raise FileNotFoundError(
                        f'GEMMA GWAS failed: {gemma_output_assoc_file}')

        # 解析结果
        gwas_df_h5 = f'{self.work_dir}/{self.name}.gemma.gwas_df.h5'
        if not have_file(gwas_df_h5):
            gwas_df = gemma_result_parser(
                gemma_output_assoc_file, f'{self.input_vcf_file}.var_stat.h5', self.genome_chr_list)
            gwas_df.to_hdf(
                gwas_df_h5, key='df', mode='w')

        self.gwas_df_h5 = gwas_df_h5


if __name__ == '__main__':
    # input vcf file
    # A VCF file with genotypes, sample order should be the same as other files
    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/landraces_indel.vcf.gz"
    # phenotype_file
    # A single column CSV file without header, values as phenotype, sample order should be the same as vcf_file
    phenotype_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/pheno.bimbam"
    # kinship_file
    # A CSV matrix without header and index, values as kinship values, sample order should be the same as phenotype_file
    kinship_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test/kinship.matrix.csv"

    # running
    # GEMMA_JOB
    name = "ai"
    work_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/test"
    mkdir(work_dir)

    gemma_job = GEMMA_JOB(name=name, work_dir=work_dir, phenotype_file=phenotype_file,
                          input_vcf_file=input_vcf_file, kinship_file=kinship_file)
    gemma_job.run()
    gemma_results_df = pd.read_hdf(gemma_job.gwas_df_h5, 'df')
