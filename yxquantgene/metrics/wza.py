"""
The WZA: A window-based method for characterizing genotype–environment associations
https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13768
"""

import numpy as np
import pandas as pd
from numpy.polynomial.polynomial import Polynomial
from scipy.stats import norm
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests
from yxutil import multiprocess_running


def calculate_empirical_p_values(stats):
    """
    Calculate empirical p-values based on genome-wide distribution.
    :param stats: numpy array of summary statistics (e.g., p-values or Bayes factors)
    :return: numpy array of empirical p-values
    """
    # Rank the statistics in ascending order
    # sorted_stats = np.sort(stats)
    ranks = np.argsort(np.argsort(stats)) + 1  # 1-based ranking
    empirical_p_values = ranks / (len(stats) + 1)  # Empirical p-value
    return empirical_p_values


def empirical_p_to_z(empirical_p_values):
    """
    Convert empirical p-values to z-scores.
    :param empirical_p_values: numpy array of empirical p-values
    :return: numpy array of z-scores
    """
    # Convert p-values to z-scores using the inverse CDF of the normal distribution
    z_scores = norm.ppf(1 - empirical_p_values)
    return z_scores


def calculate_weights(minor_allele_freq):
    """
    Calculate the weighted z-score.
    :param z_scores: numpy array of z-scores
    :param minor_allele_freq: numpy array of minor allele frequencies
    :return: weighted z-score
    """
    major_allele_freq = 1 - minor_allele_freq
    weights = minor_allele_freq * major_allele_freq
    return weights


def add_weights_to_gwas_df(gwas_df):
    """
    Add weights to the GWAS DataFrame.
    :param gwas_df: pandas DataFrame of GWAS summary statistics
        gwas_df should contain the following columns:
        - 'CHROM': chromosome IDs
        - 'POS': SNP positions
        - 'P_VALUE': p-values
        - 'MAF': minor allele frequencies
        Note: gwas_df should have all genome-wide variants.
    :return: pandas DataFrame of GWAS summary statistics with weights
    """
    gwas_df['EMPIRICAL_P_VALUE'] = calculate_empirical_p_values(
        gwas_df['P_VALUE'])
    gwas_df['Z_SCORE'] = empirical_p_to_z(gwas_df['EMPIRICAL_P_VALUE'])
    gwas_df['WEIGHT'] = calculate_weights(gwas_df['MAF'])
    return gwas_df


def get_win_z_score(gwas_df, chr_id, start_pos, end_pos):
    """
    Calculate the z-score for a given window.
    :param gwas_df: pandas DataFrame of GWAS summary statistics with weights
        gwas_df should contain the following columns:
        - 'CHROM': chromosome IDs
        - 'POS': SNP positions
        - 'WEIGHT': weights
        - 'Z_SCORE': z-scores
    """
    win_df = gwas_df[(gwas_df['CHROM'] == chr_id) & (
        gwas_df['POS'] >= start_pos) & (gwas_df['POS'] <= end_pos)]
    # win_df = gwas_df[(gwas_df['POS'] >= start_pos) & (gwas_df['POS'] <= end_pos)]
    if len(win_df) == 0:
        return None
    weighted_z_score = np.sum(
        win_df['WEIGHT'] * win_df['Z_SCORE']) / np.sqrt(np.sum(win_df['WEIGHT'] ** 2))
    return weighted_z_score, len(win_df)


def get_win_z_score_from_chr_gwas_df(chr_gwas_df, start_pos, end_pos):
    """
    Calculate the z-score for a given window.
    :param gwas_df: pandas DataFrame of GWAS summary statistics with weights
        chr_gwas_df should contain the following columns, and all SNPs should be on the same chromosome:
        - 'POS': SNP positions
        - 'WEIGHT': weights
        - 'Z_SCORE': z-scores
    """
    win_df = chr_gwas_df[(chr_gwas_df['POS'] >= start_pos)
                         & (chr_gwas_df['POS'] <= end_pos)]
    if len(win_df) == 0:
        return None
    weighted_z_score = np.sum(
        win_df['WEIGHT'] * win_df['Z_SCORE']) / np.sqrt(np.sum(win_df['WEIGHT'] ** 2))
    return weighted_z_score, len(win_df)


def get_win_z_score_for_all_gene_in_one_chr(chr_gwas_df, chr_gene_df):
    """
    Calculate the z-score for all genes in a chromosome.
    :param gwas_df: pandas DataFrame of GWAS summary statistics with weights
        gwas_df should contain the following columns:
        - 'CHROM': chromosome IDs
        - 'POS': SNP positions
        - 'WEIGHT': weights
        - 'Z_SCORE': z-scores
    :param gene_df: pandas DataFrame of gene coordinates
        gene_df should contain the following columns:
        - 'GENE': gene IDs
        - 'CHROM': chromosome IDs
        - 'START': gene start positions
        - 'END': gene end positions
    :param chr_id: chromosome ID
    :return: dictionary of gene IDs and their corresponding z-scores
    """
    z_score_list = []
    snp_num_list = []
    gene_id_list = []
    chr_id_list = []
    start_list = []
    end_list = []
    for gene_id, chr_id, start_pos, end_pos in chr_gene_df[['GENE', 'CHROM', 'START', 'END']].values:
        gene_id_list.append(gene_id)
        chr_id_list.append(chr_id)
        start_list.append(start_pos)
        end_list.append(end_pos)
        result = get_win_z_score_from_chr_gwas_df(
            chr_gwas_df, start_pos, end_pos)
        if result is not None:
            z_score, snp_num = result
            z_score_list.append(z_score)
            snp_num_list.append(snp_num)
        else:
            z_score_list.append(None)
            snp_num_list.append(0)

    output_df = pd.DataFrame()
    output_df['GENE'] = gene_id_list
    output_df['CHROM'] = chr_id_list
    output_df['START'] = start_list
    output_df['END'] = end_list
    output_df['Z_SCORE'] = z_score_list
    output_df['SNP_NUM'] = snp_num_list

    return output_df


def fit_wza_correction(gene_zw, gene_snp_count, degree=2):
    """
    拟合SNP数量与Zw的均值和标准差的多项式模型。
    :param gene_zw: 每个基因的Zw值
    :param gene_snp_count: 每个基因的SNP数量
    :param degree: 多项式拟合的阶数
    :return: 拟合的均值和标准差多项式模型
    """
    # 拟合均值模型
    mean_poly = Polynomial.fit(gene_snp_count, gene_zw, deg=degree)

    # 计算残差，并拟合标准差模型
    residuals = gene_zw - mean_poly(gene_snp_count)
    sd_poly = Polynomial.fit(gene_snp_count, np.abs(residuals), deg=degree)

    return mean_poly, sd_poly


def standardize_zw(gene_zw, gene_snp_count, mean_poly, sd_poly):
    """
    标准化Zw值。
    :param gene_zw: 每个基因的Zw值
    :param gene_snp_count: 每个基因的SNP数量
    :param mean_poly: 拟合的均值多项式模型
    :param sd_poly: 拟合的标准差多项式模型
    :return: 标准化后的Zw值
    """
    predicted_mean = mean_poly(gene_snp_count)
    predicted_sd = sd_poly(gene_snp_count)
    standardized_zw = (gene_zw - predicted_mean) / predicted_sd
    return standardized_zw


def calculate_p_values(standardized_zw):
    """
    计算p值。
    :param standardized_zw: 标准化后的Zw值
    :return: 每个基因的p值
    """
    # p_values = 2 * norm.sf(np.abs(standardized_zw))  # 双尾p值
    p_values = norm.sf(standardized_zw)  # 单尾p值
    return p_values


def adjust_zw_by_snp_num(gene_df, figsave=None):
    """
    根据SNP数量调整Zw值。
    :param gene_df: 包含基因ID、SNP数量、Zw值的DataFrame
        'GENE': 基因ID
        'SNP_NUM': SNP数量
        'Z_SCORE': Zw值
    :return: 调整后的Zw值
    """
    # 拟合SNP数量与Zw的均值和标准差的多项式模型
    mean_poly, sd_poly = fit_wza_correction(
        gene_df['Z_SCORE'], gene_df['SNP_NUM'])

    # 标准化Zw值
    standardized_zw = standardize_zw(
        gene_df['Z_SCORE'], gene_df['SNP_NUM'], mean_poly, sd_poly)
    standardized_mean_poly, standardized_sd_poly = fit_wza_correction(
        standardized_zw, gene_df['SNP_NUM'])

    # 计算p值
    p_values = calculate_p_values(standardized_zw)

    # 可视化过程
    plt.figure(figsize=(12, 6))
    gene_snp_count = gene_df['SNP_NUM']
    gene_zw = gene_df['Z_SCORE']

    # 校正前的Zw与SNP数量关系
    plt.subplot(1, 2, 1)
    plt.scatter(gene_snp_count, gene_zw, alpha=0.6, label="Original Zw")
    plt.plot(np.sort(gene_snp_count), mean_poly(
        np.sort(gene_snp_count)), color='red', label="Fitted Mean")
    plt.fill_between(
        np.sort(gene_snp_count),
        mean_poly(np.sort(gene_snp_count)) - sd_poly(np.sort(gene_snp_count)),
        mean_poly(np.sort(gene_snp_count)) + sd_poly(np.sort(gene_snp_count)),
        color='orange', alpha=0.2, label="Mean ± SD"
    )
    plt.xlabel("SNP Count")
    plt.ylabel("Zw")
    plt.title("Original Zw vs SNP Count")
    plt.legend()
    plt.grid()

    # 校正后的Zw与SNP数量关系
    plt.subplot(1, 2, 2)
    plt.scatter(gene_snp_count, standardized_zw,
                alpha=0.6, label="Standardized Zw")
    plt.plot(np.sort(gene_snp_count), standardized_mean_poly(
        np.sort(gene_snp_count)), color='red', label="Fitted Mean")
    plt.fill_between(
        np.sort(gene_snp_count),
        standardized_mean_poly(np.sort(gene_snp_count)) -
        standardized_sd_poly(np.sort(gene_snp_count)),
        standardized_mean_poly(np.sort(gene_snp_count)) +
        standardized_sd_poly(np.sort(gene_snp_count)),
        color='orange', alpha=0.2, label="Mean ± SD"
    )
    plt.xlabel("SNP Count")
    plt.ylabel("Standardized Zw")
    plt.title("Standardized Zw vs SNP Count")
    plt.legend()
    plt.grid()

    if figsave is not None:
        plt.savefig(figsave)

    plt.tight_layout()
    plt.show()
    plt.close()

    gene_df['STANDARDIZED_ZW'] = standardized_zw
    gene_df['P_VALUE'] = p_values

    return gene_df


def gene_wza_pipeline(gwas_df, gene_df, figsave=None, n_jobs=20):
    """
    WZA pipeline for gene-based analysis.
    :param gwas_df: pandas DataFrame of GWAS summary statistics
        gwas_df should contain the following columns:
        - 'CHROM': chromosome IDs
        - 'POS': SNP positions
        - 'P_VALUE': p-values
        - 'MAF': minor allele frequencies
    :param gene_df: pandas DataFrame of gene coordinates
        gene_df should contain the following columns:
        - 'GENE': gene IDs
        - 'CHROM': chromosome IDs
        - 'START': gene start positions
        - 'END': gene end positions
    :return: pandas DataFrame of gene-based analysis results
    """
    print("Adding weights to GWAS summary statistics...")
    gwas_df['EMPIRICAL_P_VALUE'] = calculate_empirical_p_values(
        gwas_df['P_VALUE'])
    gwas_df['Z_SCORE'] = empirical_p_to_z(gwas_df['EMPIRICAL_P_VALUE'])
    gwas_df['WEIGHT'] = calculate_weights(gwas_df['MAF'])
    # gwas_df.sort_values('P_VALUE')

    print("Running WZA for all genes...")
    if n_jobs == 1:
        chr_gene_df_list = []
        for chr_id in gwas_df['CHROM'].unique():
            chr_gene_df = get_win_z_score_for_all_gene_in_one_chr(gwas_df[gwas_df['CHROM'] == chr_id], gene_df[gene_df['CHROM'] == chr_id])
            chr_gene_df_list.append(chr_gene_df)
        gene_df = pd.concat(chr_gene_df_list)
        gene_df = gene_df.dropna()
    elif n_jobs > 1:
        args_dict = {}
        for chr_id in gwas_df['CHROM'].unique():
            args_dict[chr_id] = (gwas_df[gwas_df['CHROM'] == chr_id],
                                gene_df[gene_df['CHROM'] == chr_id])
        results = multiprocess_running(
            get_win_z_score_for_all_gene_in_one_chr, args_dict, n_jobs)

        gene_df = pd.concat([results[chr_id]['output']
                            for chr_id in gwas_df['CHROM'].unique()])
        gene_df = gene_df.dropna()

    print("Adjusting Zw values by SNP number...")
    gene_df = adjust_zw_by_snp_num(gene_df, figsave=figsave)

    print("Calculating FDR...")
    gene_df['FDR'] = multipletests(gene_df['P_VALUE'], method='fdr_bh')[1]

    return gene_df


if __name__ == '__main__':
    import pandas as pd
    from yxseq import read_gff_file

    gwas_df_h5 = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/gemma_output.assoc.h5'
    gff_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3'

    gwas_df = pd.read_hdf(gwas_df_h5, key='df')
    gwas_df = gwas_df.rename(columns={
                             'pvalue': 'P_VALUE', 'POS': 'POS', 'CHROM': 'CHROM', 'ID': 'ID', 'MAF': 'MAF'})

    gene_dict = read_gff_file(gff_file)['gene']
    flank_length = 2000
    gene_range_dict = {gene_id: (gene_dict[gene_id].chr_id, gene_dict[gene_id].start -
                                 flank_length, gene_dict[gene_id].end + flank_length) for gene_id in gene_dict}
    gene_id_list = list(gene_range_dict.keys())
    gene_df = pd.DataFrame()
    gene_df['GENE'] = gene_id_list
    gene_df['CHROM'] = [gene_range_dict[gene_id][0]
                         for gene_id in gene_id_list]
    gene_df['START'] = [gene_range_dict[gene_id][1]
                        for gene_id in gene_id_list]
    gene_df['END'] = [gene_range_dict[gene_id][2] for gene_id in gene_id_list]

    gene_df = gene_wza_pipeline(gwas_df, gene_df)
