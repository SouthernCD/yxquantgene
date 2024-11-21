import numpy as np
import pandas as pd
from yxmath.split import split_sequence_to_bins, bin_index_to_range, cover_bin_index
from yxutil import mkdir, log_print, multiprocess_running
from yxquantgene.utils.vcf import get_genotype_matrix_from_vcf, get_chr_list_from_var_stat_h5

"""
This module provides functions to calculate Linkage Disequilibrium (LD) matrix.
"""


def calculate_LD(query_genotype_matrix, subject_genotype_matrix, query_pos_list, subject_pos_list, ld_matrix_h5=None):
    """
    Calculate LD between each row of query_genotype_matrix and each row of subject_genotype_matrix.

    # 示例用法
    query_genotype_matrix = np.array([[1, 0, 0, 0],
                                    [0, 1, 0, 1],
                                    [0, 0, 0, 1],
                                    [0, 0, 1, 1]])
    subject_genotype_matrix = np.array([[0, 0, 1, 1],
                                        [0, 1, 0, 1],
                                        [0, 0, 0, 1],
                                        [0, 1, 1, 1],
                                        [0, 0, 1, 1]])    
    """
    # 计算相关系数矩阵
    combined_matrix = np.vstack(
        (query_genotype_matrix, subject_genotype_matrix))
    correlation_matrix = np.corrcoef(combined_matrix)

    # 提取所需的相关系数
    n_query = query_genotype_matrix.shape[0]
    n_subject = subject_genotype_matrix.shape[0]

    # 相关系数矩阵的前 n_query 行和后 n_subject 列对应的子矩阵
    result_matrix = correlation_matrix[:n_query, n_query:n_query + n_subject]

    ld_matrix = result_matrix ** 2

    if ld_matrix_h5 is not None:
        ld_df = pd.DataFrame(ld_matrix, index=query_pos_list,
                             columns=subject_pos_list)
        ld_df.to_hdf(ld_matrix_h5, key='ld_matrix')
        return ld_matrix_h5
    else:
        return ld_matrix


def get_range_genotype_matrix(genotype_matrix, var_df, start, end):
    """
    Get genotype matrix from a VCF file.
    """
    range_var_df = var_df[(var_df['POS'] >= start) & (var_df['POS'] <= end)]
    range_genotype_matrix = genotype_matrix[range_var_df.index]
    return range_genotype_matrix, list(range_var_df['POS'])


def build_LD_db(input_vcf_file, var_stat_h5_file, ld_output_dir, window_size=150000):
    mkdir(ld_output_dir, False)
    chr_list = get_chr_list_from_var_stat_h5(var_stat_h5_file)

    for chr_id in chr_list:
        log_print(f'Processing {chr_id}')
        var_df = pd.read_hdf(var_stat_h5_file, key=chr_id)
        genotype_matrix = get_genotype_matrix_from_vcf(input_vcf_file, chr_id)

        if len(var_df) == 0:
            continue

        chr_len = var_df.iloc[-1]['POS'].astype(int)

        # split sequence to bins
        w_idx_list = [i for i, s, e in split_sequence_to_bins(
            chr_len, window_size, start=1)]

        win_pair_list = []
        for w_idx in w_idx_list:
            l_w_idx = w_idx - 1
            r_w_idx = w_idx + 1
            row_idx = [w_idx]
            if l_w_idx >= 0:
                row_idx = [l_w_idx] + row_idx
            if r_w_idx <= len(w_idx_list) - 1:
                row_idx = row_idx + [r_w_idx]

            for r_idx in row_idx:
                i, j = sorted([w_idx, r_idx])
                win_pair_list.append((i, j))

        # build LD matrix
        win_pair_list = sorted(list(set(win_pair_list)))
        ld_dict = {}
        num = 0
        for q_idx, s_idx in win_pair_list:
            q_s, q_e = bin_index_to_range(q_idx, window_size, start=1)
            s_s, s_e = bin_index_to_range(s_idx, window_size, start=1)
            query_genotype_matrix, query_pos_list = get_range_genotype_matrix(
                genotype_matrix, var_df, q_s, q_e)
            subject_genotype_matrix, subject_pos_list = get_range_genotype_matrix(
                genotype_matrix, var_df, s_s, s_e)
            ld_matrix_h5 = f'{ld_output_dir}/{chr_id}_{q_idx}_{s_idx}.ld_matrix.h5'
            ld_dict[(q_idx, s_idx)] = calculate_LD(
                query_genotype_matrix, subject_genotype_matrix, query_pos_list, subject_pos_list, ld_matrix_h5)
            num += 1
            log_print(
                f'Processing {chr_id} {num}/{len(win_pair_list)} {num/len(win_pair_list) * 100:.3f}%')


def get_LD_from_db(chr_id, pos1, pos2, db_win_size, ld_db_dir):
    """
    Get LD between two positions.
    """
    q_idx = cover_bin_index(pos1, db_win_size, start=1)
    s_idx = cover_bin_index(pos2, db_win_size, start=1)

    q_idx, s_idx = sorted([q_idx, s_idx])
    ld_matrix_h5 = f'{ld_db_dir}/{chr_id}_{q_idx}_{s_idx}.ld_matrix.h5'
    ld_df = pd.read_hdf(ld_matrix_h5, key='ld_matrix')
    pos1, pos2 = sorted([pos1, pos2])
    ld = ld_df.loc[pos1, pos2]

    return ld


def get_LD_for_pairlist_from_db(chr_id, pos_pair_list, db_win_size, ld_db_dir):
    """
    Get LD for a list of position pairs.
    """
    pos_pair_list = [(pos1, pos2) if pos1 < pos2 else (pos2, pos1)
                     for pos1, pos2 in pos_pair_list]
    # pos_pair_dict = {(pos1, pos2): None if pos1 < pos2 else (
    #     pos2, pos1) for pos1, pos2 in pos_pair_list}
    pos_pair_dict = {(pos1, pos2): (cover_bin_index(pos1, db_win_size, start=1), cover_bin_index(
        pos2, db_win_size, start=1))for pos1, pos2 in pos_pair_list}
    ld_idx_list = list(set(pos_pair_dict.values()))

    ld_df_dict = {}
    for ld_id_pair in ld_idx_list:
        ld_df_dict[ld_id_pair] = pd.read_hdf(
            f'{ld_db_dir}/{chr_id}_{ld_id_pair[0]}_{ld_id_pair[1]}.ld_matrix.h5', key='ld_matrix')

    ld_dict = {p: ld_df_dict[(
        pos_pair_dict[p][0], pos_pair_dict[p][1])].loc[p[0], p[1]] for p in pos_pair_dict}

    return ld_dict


def get_LD_decay_for_one_win(win_idx, flank_win_idx_list, var_pos_idx_df, ld_db_path, chr_id, max_decay_size):
    """
    windows size of LD database have to bigger than half of max_decay_size
    var_pos_idx_df = var_df.reset_index().set_index('POS')
    var_pos_idx_df
    ld_db_path = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/snp_ld"
    chr_id = "Chr01"
    win_idx = 1
    flank_win_idx_list = [1, 2]
    max_decay_size = 500000        
    """

    # 读取左中右三个窗口的 LD 矩阵
    ld_df_list = []
    for r_idx in flank_win_idx_list:
        if r_idx >= win_idx:
            ld_chunk_matrix_h5 = f'{ld_db_path}/{chr_id}_{win_idx}_{r_idx}.ld_matrix.h5'
            ld_df = pd.read_hdf(ld_chunk_matrix_h5, key='ld_matrix')
            ld_df = ld_df.T
        else:
            ld_chunk_matrix_h5 = f'{ld_db_path}/{chr_id}_{r_idx}_{win_idx}.ld_matrix.h5'
            ld_df = pd.read_hdf(ld_chunk_matrix_h5, key='ld_matrix')
            ld_df
        ld_df_list.append(ld_df)
    # 合并 LD 矩阵，并筛选出有效行和列
    ld_df = pd.concat(ld_df_list, axis=0)
    ld_df = ld_df.T

    valid_rows = ld_df.index.intersection(var_pos_idx_df.index)
    valid_cols = ld_df.columns.intersection(var_pos_idx_df.index)
    ld_df = ld_df.loc[valid_rows, valid_cols]

    # 遍历 LD 矩阵的每一行（目标窗口中的所有变异位点）

    dist_vs_ld_df_list = []
    num = 0
    for q_pos in ld_df.index:
        q_pos_ld_df = pd.DataFrame(
            {'dist': ld_df.loc[q_pos].index - q_pos, 'ld': ld_df.loc[q_pos].values})
        q_pos_ld_df = q_pos_ld_df[(q_pos_ld_df['dist'] > 0) & (
            q_pos_ld_df['dist'] < max_decay_size)]
        dist_vs_ld_df_list.append(q_pos_ld_df)
        num += 1
        # print(f"Processing {chr_id} {win_idx} {num}/{len(ld_df.index)} {num/len(ld_df.index) * 100:.3f}%")

    dist_vs_ld_df = pd.concat(dist_vs_ld_df_list, ignore_index=True)

    dist_vs_ld_df[(dist_vs_ld_df['dist'] <= 450) & (
        dist_vs_ld_df['dist'] >= 350)]['ld'].median()

    return dist_vs_ld_df


def get_chr_LD_dist_df(chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size=150000, max_decay_size=150000, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, threads=20):
    # 读取变异位点信息
    var_df = pd.read_hdf(var_stat_h5_file, key=chr_id)
    # prune variants based on var_stat
    if max_missing_rate is not None:
        var_df = var_df[(var_df['MISSF'] <= max_missing_rate)]
    if min_maf is not None:
        var_df = var_df[(var_df['MAF'] >= min_maf)]
    if max_het_rate is not None:
        var_df = var_df[(var_df['HETF'] <= max_het_rate)]

    # 获取染色体长度
    chr_len = var_df.iloc[-1]['POS'].astype(int)

    # 将变异位点信息重置索引并设置 POS 为索引
    var_pos_idx_df = var_df.reset_index().set_index('POS')

    # 将染色体序列分割成窗口
    w_idx_list = [i for i, s, e in split_sequence_to_bins(
        chr_len, ld_db_win_size, start=1)]

    # 遍历每个窗口
    args_dict = {}
    dist_vs_ld_df = pd.DataFrame({'dist': [], 'ld': []})
    for w_idx in w_idx_list:
        # 获取右窗口索引
        # l_w_idx = w_idx - 1
        r_w_idx = w_idx + 1
        row_idx = [w_idx]
        # if l_w_idx >= 0:
        #     row_idx = [l_w_idx] + row_idx
        if r_w_idx <= len(w_idx_list) - 1:
            row_idx = row_idx + [r_w_idx]
        args_dict[w_idx] = (w_idx, row_idx, var_pos_idx_df,
                            ld_db_path, chr_id, max_decay_size)
        # log_print(
        #     f"Processing {chr_id} {w_idx}/{len(w_idx_list)} {w_idx/len(w_idx_list) * 100:.3f}%")

        # dist_vs_ld_df = pd.concat([dist_vs_ld_df, get_LD_decay_for_one_win(
        #     w_idx, row_idx, var_pos_idx_df, ld_db_path, chr_id, max_decay_size)], ignore_index=True)

    # 并行运行
    mlt_dict = multiprocess_running(
        get_LD_decay_for_one_win, args_dict, threads)

    dist_vs_ld_df_list = [mlt_dict[i]['output'] for i in mlt_dict]
    dist_vs_ld_df = pd.concat(dist_vs_ld_df_list, ignore_index=True)

    return dist_vs_ld_df


def get_chr_LD_decay_curve(chr_id, var_stat_h5_file, ld_db_path, bin_size=100, ld_db_win_size=150000, max_decay_size=150000, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, threads=20):
    dist_vs_ld_df = get_chr_LD_dist_df(chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size,
                                       max_decay_size, max_missing_rate, min_maf, max_het_rate, threads)

    ld_decay_df = pd.DataFrame({'start': [], 'end': [], 'mean': [], 'median': [], 'std': []})
    for i, s, e in split_sequence_to_bins(max_decay_size, bin_size, start=1):
        print(f"Processing {s}-{e}")
        bin_dist_vs_ld_df = dist_vs_ld_df[(dist_vs_ld_df['dist'] <= e) & (
            dist_vs_ld_df['dist'] >= s)]
        median = bin_dist_vs_ld_df['ld'].median()
        mean = bin_dist_vs_ld_df['ld'].mean()
        std = bin_dist_vs_ld_df['ld'].std()
        ld_decay_df = pd.concat([ld_decay_df, pd.DataFrame({'start': [s], 'end': [e], 'mean': [mean], 'median': [median], 'std': [std]})])

    return ld_decay_df

if __name__ == '__main__':

    chr_id = "Chr01"
    var_stat_h5_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/landraces_snp_stat.h5"
    ld_db_path = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/snp_ld"
    ld_db_win_size = 500000
    max_decay_size = 500000
    bin_size = 100
    max_missing_rate = 0.5
    min_maf = 0.01
    max_het_rate = 0.1
    threads = 20
