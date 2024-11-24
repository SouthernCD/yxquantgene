import numpy as np
import pandas as pd
from yxmath.split import split_sequence_to_bins, bin_index_to_range, cover_bin_index
from yxutil import mkdir, log_print, multiprocess_running
from yxquantgene.utils.vcf import get_genotype_matrix_from_vcf, get_chr_list_from_var_stat_h5
from scipy import interpolate

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


def get_ld_df_sum(df1, df2):
    merged_df = pd.merge(df1, df2, on='bin', suffixes=('_df1', '_df2'))
    # 对应 bin 的 ld_sum 和 count 进行加和
    merged_df['ld_sum_total'] = merged_df['ld_sum_df1'] + \
        merged_df['ld_sum_df2']
    merged_df['count_total'] = merged_df['count_df1'] + merged_df['count_df2']

    # 选择需要的列
    result_df = merged_df[['bin', 'ld_sum_total', 'count_total']]
    result_df = result_df.rename(
        columns={'ld_sum_total': 'ld_sum', 'count_total': 'count'})

    return result_df


def get_ld_df_list_sum(ld_df_list):
    ld_df_sum = ld_df_list[0]
    for ld_df in ld_df_list[1:]:
        ld_df_sum = get_ld_df_sum(ld_df_sum, ld_df)
    return ld_df_sum


def get_LD_decay_mean_for_one_win(win_idx, flank_win_idx_list, var_pos_idx_df, ld_db_path, chr_id, max_decay_size, curve_bin_size):
    """
    windows size of LD database have to bigger than half of max_decay_size
    var_pos_idx_df = var_df.reset_index().set_index('POS')
    var_pos_idx_df
    ld_db_path = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/snp_ld"
    chr_id = "Chr01"
    win_idx = 1
    flank_win_idx_list = [1, 2]
    max_decay_size = 500000
    curve_bin_size = 100        
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

    bin_stats = get_LD_decay_mean_from_ld_df(ld_df, max_decay_size=max_decay_size, curve_bin_size=curve_bin_size)

    return bin_stats


def get_LD_decay_mean_from_ld_df(ld_df, max_decay_size=None, curve_bin_size=100):
    """
    ld_df: DataFrame with index and columns as variant positions, and values as LD values.
    max_decay_size: 最大 LD 衰减距离。
    curve_bin_size: 曲线的 bin 大小。
    """
    if max_decay_size is None:
        max_decay_size = ld_df.index.max() - ld_df.index.min() + 1

    # 遍历 LD 矩阵的每一行（目标窗口中的所有变异位点）
    bins = pd.interval_range(start=0, end=max_decay_size,
                             freq=curve_bin_size, closed='right')
    bin_stats = pd.DataFrame({
        'bin': bins,   # 使用区间作为 bin 列
        'ld_sum': 0,   # ld_sum 列初始化为 0
        'count': 0     # count 列初始化为 0
    })

    num = 0
    q_pos_ld_df_list = []
    for q_pos in ld_df.index:
        q_pos_ld_df = pd.DataFrame(
            {'dist': ld_df.loc[q_pos].index - q_pos, 'ld': ld_df.loc[q_pos].values})
        q_pos_ld_df = q_pos_ld_df[(q_pos_ld_df['dist'] > 0) & (
            q_pos_ld_df['dist'] < max_decay_size)]
        q_pos_ld_df_list.append(q_pos_ld_df)
        num += 1

        if len(q_pos_ld_df_list) > 1000:
            # break
            q_pos_ld_df = pd.concat(q_pos_ld_df_list, ignore_index=True)
            q_pos_ld_df = q_pos_ld_df.sort_values('dist')
            bins = range(0, max_decay_size + curve_bin_size, curve_bin_size)
            q_pos_ld_df['bin'] = pd.cut(
                q_pos_ld_df['dist'], bins=bins, right=True)
            bin_stats_now = q_pos_ld_df.groupby('bin').agg(
                ld_sum=('ld', 'sum'),  # 计算 ld 的总和
                count=('ld', 'count')  # 计算 count
            ).reset_index()
            bin_stats = get_ld_df_sum(bin_stats, bin_stats_now)
            q_pos_ld_df_list = []

            # print(f"Processing {chr_id} {win_idx} {num}/{len(ld_df.index)} {num/len(ld_df.index) * 100:.3f}%")

    if len(q_pos_ld_df_list) > 0:
        q_pos_ld_df = pd.concat(q_pos_ld_df_list, ignore_index=True)
        q_pos_ld_df = q_pos_ld_df.sort_values('dist')
        bins = range(0, max_decay_size + curve_bin_size, curve_bin_size)
        q_pos_ld_df['bin'] = pd.cut(q_pos_ld_df['dist'], bins=bins, right=True)
        bin_stats_now = q_pos_ld_df.groupby('bin').agg(
            ld_sum=('ld', 'sum'),  # 计算 ld 的总和
            count=('ld', 'count')  # 计算 count
        ).reset_index()
        bin_stats = get_ld_df_sum(bin_stats, bin_stats_now)
        q_pos_ld_df_list = []

        # print(f"Processing {chr_id} {win_idx} {num}/{len(ld_df.index)} {num/len(ld_df.index) * 100:.3f}%")

    return bin_stats

def get_LD_decay_mean_for_one_chr(chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size=150000, max_decay_size=150000, curve_bin_size=100, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, threads=20):
    # 读取变异位点信息
    var_df = pd.read_hdf(var_stat_h5_file, key=chr_id)

    if len(var_df) == 0:
        return None

    # prune variants based on var_stat
    if max_missing_rate is not None:
        var_df = var_df[(var_df['MISSF'] <= max_missing_rate)]
    if min_maf is not None:
        var_df = var_df[(var_df['MAF'] >= min_maf)]
    if max_het_rate is not None:
        var_df = var_df[(var_df['HETF'] <= max_het_rate)]

    if len(var_df) == 0:
        return None

    # 获取染色体长度
    chr_len = var_df.iloc[-1]['POS'].astype(int)

    # 将变异位点信息重置索引并设置 POS 为索引
    var_pos_idx_df = var_df.reset_index().set_index('POS')

    # 将染色体序列分割成窗口
    w_idx_list = [i for i, s, e in split_sequence_to_bins(
        chr_len, ld_db_win_size, start=1)]

    # 遍历每个窗口
    args_dict = {}
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
                            ld_db_path, chr_id, max_decay_size, curve_bin_size)

        # log_print(
        #     f"Processing {chr_id} {w_idx}/{len(w_idx_list)} {w_idx/len(w_idx_list) * 100:.3f}%")

        # dist_vs_ld_df = pd.concat([dist_vs_ld_df, get_LD_decay_for_one_win(
        #     w_idx, row_idx, var_pos_idx_df, ld_db_path, chr_id, max_decay_size)], ignore_index=True)

    # 并行运行
    mlt_dict = multiprocess_running(
        get_LD_decay_mean_for_one_win, args_dict, threads)

    # 合并结果
    ld_df_list = [mlt_dict[i]['output'] for i in mlt_dict]
    ld_df_sum = get_ld_df_list_sum(ld_df_list)

    ld_df_sum['ld_mean'] = ld_df_sum['ld_sum'] / ld_df_sum['count']

    ld_decay_df = ld_df_sum[['bin', 'ld_mean']]

    return ld_decay_df


def get_LD_decay(var_stat_h5_file, ld_db_path, ld_db_win_size=150000, max_decay_size=150000, curve_bin_size=100, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, threads=20):

    chr_list = get_chr_list_from_var_stat_h5(var_stat_h5_file)
    ld_decay_df_dict = {}
    for chr_id in chr_list:
        print(f'Processing {chr_id}')
        ld_decay_df = get_LD_decay_mean_for_one_chr(
            chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size, max_decay_size, curve_bin_size, max_missing_rate, min_maf, max_het_rate, threads)
        if ld_decay_df is not None:
            ld_decay_df_dict[chr_id] = ld_decay_df

    # 把chr_id添加为一列，并整合到一个DataFrame中
    ld_decay_df = pd.concat(ld_decay_df_dict, keys=ld_decay_df_dict.keys(), names=[
                            'chr_id']).reset_index(level=0)

    return ld_decay_df


def get_half_ld_dist(ld_decay_df, half_ld_value=0.5):
    ld_decay_df['bin'] = ld_decay_df['bin'].apply(lambda x: x.right).values
    
    if 'ld_mean' not in ld_decay_df.columns:
        ld_decay_df['ld_mean'] = ld_decay_df['ld_sum'] / ld_decay_df['count']

    
    ld_decay_df = ld_decay_df.dropna()

    ld_decay_df = ld_decay_df.sort_values('bin')

    if all(ld_decay_df['ld_mean'] < half_ld_value):
        return ld_decay_df['bin'].min()
    elif all(ld_decay_df['ld_mean'] > half_ld_value):
        return ld_decay_df['bin'].max()
    else:
        ld2dist_func = interpolate.interp1d(ld_decay_df['ld_mean'], ld_decay_df['bin'])
        half_ld_dist = ld2dist_func(half_ld_value)
        return half_ld_dist

def get_double_ld_df(s_win_idx,e_win_idx,ld_db_path,chr_id):
    if s_win_idx == e_win_idx:
        ld_chunk_matrix_h5 = f'{ld_db_path}/{chr_id}_{s_win_idx}_{s_win_idx}.ld_matrix.h5'
        ld_df = pd.read_hdf(ld_chunk_matrix_h5, key='ld_matrix')
    else:
        ss_matrix_h5 = f'{ld_db_path}/{chr_id}_{s_win_idx}_{s_win_idx}.ld_matrix.h5'
        ss_df = pd.read_hdf(ss_matrix_h5, key='ld_matrix')
        ee_matrix_h5 = f'{ld_db_path}/{chr_id}_{e_win_idx}_{e_win_idx}.ld_matrix.h5'
        ee_df = pd.read_hdf(ee_matrix_h5, key='ld_matrix')
        se_matrix_h5 = f'{ld_db_path}/{chr_id}_{s_win_idx}_{e_win_idx}.ld_matrix.h5'
        se_df = pd.read_hdf(se_matrix_h5, key='ld_matrix')
        a_df = pd.concat([ss_df, se_df], axis=1)
        b_df = pd.concat([se_df.T, ee_df], axis=1)
        ld_df = pd.concat([a_df, b_df], axis=0)
    return ld_df

def get_half_ld_dist_for_one_chr(chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size=150000, stat_win_size=150000, max_decay_size=150000, curve_bin_size=100, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, half_ld_value=0.5):
    # 读取变异位点信息
    var_df = pd.read_hdf(var_stat_h5_file, key=chr_id)

    if len(var_df) == 0:
        return None

    # prune variants based on var_stat
    if max_missing_rate is not None:
        var_df = var_df[(var_df['MISSF'] <= max_missing_rate)]
    if min_maf is not None:
        var_df = var_df[(var_df['MAF'] >= min_maf)]
    if max_het_rate is not None:
        var_df = var_df[(var_df['HETF'] <= max_het_rate)]

    if len(var_df) == 0:
        return None

    # 获取染色体长度
    chr_len = var_df.iloc[-1]['POS'].astype(int)

    # 将变异位点信息重置索引并设置 POS 为索引
    var_pos_idx_df = var_df.reset_index().set_index('POS')

    # 将染色体序列分割成窗口
    stat_w_idx_list = [(i,s,e) for i, s, e in split_sequence_to_bins(
        chr_len, stat_win_size, start=1)]

    
    job_group_dict = {}
    for stat_w_idx, stat_s, stat_e in stat_w_idx_list:
        ld_db_s_w_idx = cover_bin_index(stat_s, ld_db_win_size, start=1)
        ld_db_e_w_idx = cover_bin_index(stat_e, ld_db_win_size, start=1)
        job_group_dict.setdefault((ld_db_s_w_idx, ld_db_e_w_idx), []).append((stat_w_idx, stat_s, stat_e))

    half_ld_dist_df = pd.DataFrame(columns=['chr_id', 'bin', 'half_ld_dist'])

    for ld_db_w_idxs in job_group_dict:
        s_win_idx, e_win_idx = ld_db_w_idxs
        # print(f"Processing {chr_id} {s_win_idx*ld_db_win_size}/{chr_len} {s_win_idx*ld_db_win_size/chr_len * 100:.3f}%")
        ld_df = get_double_ld_df(s_win_idx,e_win_idx,ld_db_path,chr_id)
        
        for stat_w_idx, stat_s, stat_e in job_group_dict[ld_db_w_idxs]:
            win_var_pos_idx_df = var_pos_idx_df[(var_pos_idx_df.index >= stat_s) & (var_pos_idx_df.index <= stat_e)]
            valid_rows = ld_df.index.intersection(win_var_pos_idx_df.index)
            valid_cols = ld_df.columns.intersection(win_var_pos_idx_df.index)
            stat_win_ld_df = ld_df.loc[valid_rows, valid_cols]
            ld_decay_df = get_LD_decay_mean_from_ld_df(stat_win_ld_df, max_decay_size=max_decay_size, curve_bin_size=curve_bin_size)
            ld_decay_df['ld_mean'] = ld_decay_df['ld_sum'] / ld_decay_df['count']
            half_ld_dist = get_half_ld_dist(ld_decay_df, half_ld_value)
            half_ld_dist_df = pd.concat([half_ld_dist_df, pd.DataFrame({'chr_id': [chr_id], 'bin': [(stat_s, stat_e)], 'half_ld_dist': [half_ld_dist]})], ignore_index=True)
            
    return half_ld_dist_df


def get_half_ld_dist_for_genome(var_stat_h5_file, ld_db_path, ld_db_win_size=150000, stat_win_size=150000, max_decay_size=150000, curve_bin_size=100, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, half_ld_value=0.5, threads=20):
    chr_list = get_chr_list_from_var_stat_h5(var_stat_h5_file)

    args_dict = {}
    for chr_id in chr_list:
        args_dict[chr_id] = (chr_id, var_stat_h5_file, ld_db_path, ld_db_win_size, stat_win_size, max_decay_size, curve_bin_size, max_missing_rate, min_maf, max_het_rate, half_ld_value)
    
    mlt_dict = multiprocess_running(get_half_ld_dist_for_one_chr, args_dict, threads)
    half_ld_dist_df_dict = {i:mlt_dict[i]['output'] for i in mlt_dict if mlt_dict[i]['output'] is not None}

    half_ld_dist_df = pd.concat(half_ld_dist_df_dict)

    return half_ld_dist_df
    


if __name__ == '__main__':

    var_stat_h5_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/landraces_snp_stat.h5"
    ld_db_path = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.reseq_GWAS/population_structure/snp_ld"
    ld_db_win_size = 500000
    max_decay_size = 150000
    stat_win_size = 100000
    curve_bin_size = 50
    threads = 20
    half_ld_value = 0.5
    max_missing_rate=0.5
    min_maf=0.1
    max_het_rate=0.01

    half_ld_dist_df = get_half_ld_dist_for_genome(var_stat_h5_file, ld_db_path, ld_db_win_size, stat_win_size, max_decay_size, curve_bin_size, max_missing_rate, min_maf, max_het_rate, half_ld_value, threads)
