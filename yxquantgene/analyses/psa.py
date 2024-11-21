from yxutil import multiprocess_running
import numpy as np
import pandas as pd
from yxquantgene.metrics.ld import get_LD_for_pairlist_from_db
from yxquantgene.utils.vcf import extract_subvcf_by_varIDs, get_chr_list_from_var_stat_h5
from yxmath.split import split_sequence_to_bins, bin_index_to_range, cover_bin_index

# def get_representative_positions_from_chunk(ld_chunk_matrix_h5, var_pos_idx_df, ld_r2_threshold):
#     # 从 HDF5 文件中读取 LD 矩阵，并筛选出有效行和列
#     ld_df = pd.read_hdf(ld_chunk_matrix_h5, key='ld_matrix')
#     valid_rows = ld_df.index.intersection(var_pos_idx_df.index)
#     valid_cols = ld_df.columns.intersection(var_pos_idx_df.index)
#     ld_df = ld_df.loc[valid_rows, valid_cols]

#     # 初始化代表性位置列表
#     rep_pos_list = []

#     # 遍历 LD 矩阵的每一行
#     for q_pos in ld_df.index:
#         # 筛选出其他位置
#         s_pos_list = [i for i in ld_df.columns if i != q_pos]

#         # 计算 LD 和 MAF 条件
#         ## 判断该chunk中的其他位置是否有比当前位置更大的 LD 和 MAF
#         s_pos_bigger_ld = ld_df.loc[q_pos][s_pos_list].values > ld_r2_threshold
#         s_pos_bigger_maf = var_pos_idx_df.loc[s_pos_list].MAF > var_pos_idx_df.loc[q_pos].MAF
#         s_pos_bigger = s_pos_bigger_ld & s_pos_bigger_maf

#         # 判断是否所有条件都不满足
#         if not s_pos_bigger.any():
#             # 添加代表性位置
#             rep_pos_list.append(q_pos)

#     # 返回代表性位置列表
#     return rep_pos_list


# def preprune_by_ld_chunk(var_df, chr_id, ld_db_path, ld_r2_threshold, window_size, threads=20):
#     """
#     根据 LD 矩阵对变异位点进行预修剪。

#     参数:
#     var_df (pd.DataFrame): 包含变异位点信息的 DataFrame。
#     chr_id (str): 染色体 ID。
#     ld_db_path (str): LD 矩阵数据库路径。
#     ld_r2_threshold (float): LD r^2 阈值。
#     window_size (int): 窗口大小。

#     返回:
#     pd.DataFrame: 修剪后的变异位点 DataFrame。
#     """
#     # 获取染色体长度
#     chr_len = var_df.iloc[-1]['POS'].astype(int)

#     # 将变异位点信息重置索引并设置 POS 为索引
#     var_pos_idx_df = var_df.reset_index().set_index('POS')

#     # 将染色体序列分割成窗口
#     w_idx_list = [i for i, s, e in split_sequence_to_bins(chr_len, window_size, start=1)]

#     # 生成窗口对列表
#     win_pair_list = []
#     for w_idx in w_idx_list:
#         l_w_idx = w_idx - 1
#         r_w_idx = w_idx + 1
#         row_idx = [w_idx]
#         if l_w_idx >= 0:
#             row_idx = [l_w_idx] + row_idx
#         if r_w_idx <= len(w_idx_list) - 1:
#             row_idx = row_idx + [r_w_idx]

#         for r_idx in row_idx:
#             i, j = sorted([w_idx, r_idx])
#             win_pair_list.append((i, j))

#     print(f'Processing {chr_id} from chunks...')

#     # 创建参数字典
#     args_dict = {}
#     for q_idx, s_idx in win_pair_list:
#         ld_chunk_matrix_h5 = f'{ld_db_path}/{chr_id}_{q_idx}_{s_idx}.ld_matrix.h5'
#         args_dict[(q_idx, s_idx)] = (ld_chunk_matrix_h5, var_pos_idx_df, ld_r2_threshold)

#     # 并行运行 get_representative_positions_from_chunk 函数
#     mlt_dict = multiprocess_running(get_representative_positions_from_chunk, args_dict, threads)

#     # 收集代表性位置
#     rep_pos_list = []
#     for q_idx, s_idx in win_pair_list:
#         rep_pos_list.extend(mlt_dict[(q_idx, s_idx)]['output'])
#     rep_pos_list = list(set(rep_pos_list))
#     rep_pos_list = sorted(rep_pos_list)

#     # 筛选修剪后的变异位点
#     rep_var_df = var_df[var_df['POS'].isin(rep_pos_list)]

#     return rep_var_df


def pruner_for_one_win(win_idx, flank_win_idx_list, ld_db_path, var_pos_idx_df, chr_id, ld_decay_size, ld_r2_threshold):
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
    rep_pos_list = []
    for q_pos in ld_df.index:
        # 筛选出其他位置
        s_pos_list = [i for i in ld_df.columns if i != q_pos and abs(i - q_pos) < ld_decay_size]

        # 计算 LD 和 MAF 条件
        ## 判断该chunk中的其他位置是否有比当前位置更大的 LD 和 MAF
        s_pos_bigger_ld = ld_df.loc[q_pos][s_pos_list].values > ld_r2_threshold
        s_pos_bigger_maf = var_pos_idx_df.loc[s_pos_list].MAF > var_pos_idx_df.loc[q_pos].MAF
        s_pos_bigger = s_pos_bigger_ld & s_pos_bigger_maf

        # 判断是否所有条件都不满足
        if not s_pos_bigger.any():
            # 添加代表性位置
            rep_pos_list.append(q_pos)

    return rep_pos_list


def pruner_by_traveling_window(chr_id, var_df, ld_db_path, ld_db_win_size=150000, ld_decay_size=150000, ld_r2_threshold=0.5, threads=20):
    # 获取染色体长度
    chr_len = var_df.iloc[-1]['POS'].astype(int)

    # 将变异位点信息重置索引并设置 POS 为索引
    var_pos_idx_df = var_df.reset_index().set_index('POS')

    # 将染色体序列分割成窗口
    w_idx_list = [i for i, s, e in split_sequence_to_bins(chr_len, ld_db_win_size, start=1)]

    # 遍历每个窗口
    args_dict = {}
    for w_idx in w_idx_list:
        # 获取左右窗口索引
        l_w_idx = w_idx - 1
        r_w_idx = w_idx + 1
        row_idx = [w_idx]
        if l_w_idx >= 0:
            row_idx = [l_w_idx] + row_idx
        if r_w_idx <= len(w_idx_list) - 1:
            row_idx = row_idx + [r_w_idx]
        args_dict[w_idx] = (w_idx, row_idx, ld_db_path, var_pos_idx_df, chr_id, ld_decay_size, ld_r2_threshold)

    # 并行运行 pruner_for_one_win 函数
    mlt_dict = multiprocess_running(pruner_for_one_win, args_dict, threads)
    rep_pos_list = []
    for w_idx in w_idx_list:
        rep_pos_list.extend(mlt_dict[w_idx]['output'])
    rep_pos_list = list(set(rep_pos_list))
    rep_pos_list = sorted(rep_pos_list)

    # 筛选修剪后的变异位点
    rep_var_df = var_df[var_df['POS'].isin(rep_pos_list)]

    return rep_var_df



def psa_snp_pruner(input_vcf_file, var_stat_h5_file, ld_db_path, output_prefix, ld_db_win_size=150000, ld_decay_size=150000, ld_r2_threshold=0.5, max_missing_rate=0.5, min_maf=0.01, max_het_rate=0.5, threads=20):
    """
    Prune variants based on LD.
    """
    output_prefix = output_prefix + f'.win{window_size}.'
    if max_missing_rate is not None:
        output_prefix = output_prefix + f'miss{max_missing_rate}.'
    if min_maf is not None:
        output_prefix = output_prefix + f'maf{min_maf}.'
    if max_het_rate is not None:
        output_prefix = output_prefix + f'het{max_het_rate}.'
    if ld_r2_threshold is not None:
        output_prefix = output_prefix + f'rq{ld_r2_threshold}.'

    chr_list = get_chr_list_from_var_stat_h5(var_stat_h5_file)

    chr_rep_var_df_dict = {}
    for chr_id in chr_list:
        # read var_stat
        var_df = pd.read_hdf(var_stat_h5_file, key=chr_id)
        if len(var_df) == 0:
            continue        

        print(f'Processing {chr_id}, {len(var_df)} variants...')

        # prune variants based on var_stat
        if max_missing_rate is not None:
            var_df = var_df[(var_df['MISSF'] <= max_missing_rate)]
        if min_maf is not None:
            var_df = var_df[(var_df['MAF'] >= min_maf)]
        if max_het_rate is not None:
            var_df = var_df[(var_df['HETF'] <= max_het_rate)]

        print(f'After QC, {len(var_df)} variants left.')

        # # prune variants based on LD chunk
        # var_df = preprune_by_ld_chunk(var_df, chr_id, ld_db_path, ld_r2_threshold, window_size, threads)

        # prune variants based on LD traveling window
        rep_var_df = pruner_by_traveling_window(chr_id, var_df, ld_db_path, ld_db_win_size=ld_db_win_size, ld_decay_size=ld_decay_size, ld_r2_threshold=ld_r2_threshold, threads=threads)

        print(f'After LD pruning, {len(rep_var_df)} variants left.')

        chr_rep_var_df_dict[chr_id] = rep_var_df

    output_var_stat_h5_file = output_prefix + '.var_stat.h5'

    for chr_id in chr_rep_var_df_dict:
        chr_rep_var_df_dict[chr_id].to_hdf(
            output_var_stat_h5_file, key=chr_id, mode='a')

    pruned_var_list = []
    for chr_id in chr_rep_var_df_dict:
        pruned_var_list.extend(chr_rep_var_df_dict[chr_id]['POS'].tolist())

    print(f'After all pruning, {len(pruned_var_list)} variants left.')

    print('Extracting sub VCF file...')
    output_vcf_file = output_prefix + '.vcf.gz'
    extract_subvcf_by_varIDs(input_vcf_file, pruned_var_list, output_vcf_file)

if __name__ == '__main__':
    from yxquantgene import build_LD_db, build_var_stat_table, psa_snp_pruner

    # Load the genotype matrix
    test_vcf_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples.vcf.gz'
    ref_genome_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa'
    var_stat_h5_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples.var_stat.h5'

    build_var_stat_table(test_vcf_file, ref_genome_file, var_stat_h5_file)

    # build LD database
    window_size = 150000

    ld_db_path = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples_ld'
    build_LD_db(test_vcf_file, var_stat_h5_file,
                ld_db_path, window_size=window_size)

    # Prune variants for population structure analysis
    window_size = 150000
    ld_r2_threshold = 0.5
    max_missing_rate = 0.5
    min_maf = 0.01
    max_het_rate = 0.5
    threads = 20

    output_prefix = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples_pruned'
    psa_snp_pruner(test_vcf_file, var_stat_h5_file, ld_db_path, output_prefix, window_size=window_size,
                    ld_r2_threshold=ld_r2_threshold, max_missing_rate=max_missing_rate, min_maf=min_maf, max_het_rate=max_het_rate, threads=threads)
