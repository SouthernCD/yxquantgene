import pandas as pd
from interlap import InterLap
import numpy as np
from yxutil import multiprocess_running


def in_region(x, f_tag, foreground_interlap_dict, b_tag, background_interlap_dict):
    chr_id = x['CHROM']
    pos = x['POS']
    if chr_id in foreground_interlap_dict and len(list(foreground_interlap_dict[chr_id].find((pos, pos)))) > 0:
        return f_tag
    elif chr_id in background_interlap_dict and len(list(background_interlap_dict[chr_id].find((pos, pos)))) > 0:
        return b_tag
    else:
        return np.nan


def add_site_tag_by_region(gwas_asso_df, tag, foreground_tag, foreground_region_dict, background_tag, background_region_dict):

    foreground_interlap_dict = {}
    for chr_id in foreground_region_dict:
        foreground_interlap_dict[chr_id] = InterLap()
        foreground_interlap_dict[chr_id].update(
            [(start, end) for start, end in foreground_region_dict[chr_id]])

    background_interlap_dict = {}
    for chr_id in background_region_dict:
        background_interlap_dict[chr_id] = InterLap()
        background_interlap_dict[chr_id].update(
            [(start, end) for start, end in background_region_dict[chr_id]])

    gwas_asso_df[tag] = gwas_asso_df.apply(lambda x: in_region(
        x, foreground_tag, foreground_interlap_dict, background_tag, background_interlap_dict), axis=1)

    return gwas_asso_df


def read_region_from_bed_file(input_bed_file):
    region_dict = {}
    with open(input_bed_file) as f:
        for l in f:
            l = l.strip().split("\t")
            if l[0] not in region_dict:
                region_dict[l[0]] = []
            region_dict[l[0]].append((int(l[1]), int(l[2])))
    return region_dict


def perm_job(input_df, test_tag, foreground_tag, background_tag, test_value):
    perm_df = input_df.copy()
    perm_df['PERM_TAG'] = np.random.permutation(perm_df[test_tag].values)
    perm_foreground = np.percentile(
        perm_df[perm_df['PERM_TAG'] == foreground_tag][test_value], 95)
    perm_background = np.percentile(
        perm_df[perm_df['PERM_TAG'] == background_tag][test_value], 95)
    perm_ratio = perm_foreground / perm_background
    return perm_ratio


def site_region_enrichment_test(gwas_asso_df, test_tag, foreground_tag, background_tag, test_value, permutation_num, threads=10):
    observe_foreground = np.percentile(
        gwas_asso_df[gwas_asso_df[test_tag] == foreground_tag][test_value], 95)
    observe_background = np.percentile(
        gwas_asso_df[gwas_asso_df[test_tag] == background_tag][test_value], 95)
    observe_ratio = observe_foreground / observe_background

    args_dict = {}
    tmp_df = gwas_asso_df[[test_tag, test_value]].copy()
    for i in range(permutation_num):
        args_dict[i] = (tmp_df, test_tag, foreground_tag,
                        background_tag, test_value)

    mlt_dict = multiprocess_running(perm_job, args_dict, threads)

    perm_ratio_list = np.array([mlt_dict[i]['output']
                               for i in range(permutation_num)])

    p_value = (np.sum(perm_ratio_list >= observe_ratio) + 1) / \
        (permutation_num + 1)

    return observe_ratio, perm_ratio_list, p_value


if __name__ == '__main__':

    open_chromatin_v5_bed_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Data/atacseq/open_chromatin.v5.bed"
    close_chromatin_v5_bed_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Data/atacseq/close_chromatin.v5.bed"

    open_chromatin_dict = read_region_from_bed_file(open_chromatin_v5_bed_file)
    close_chromatin_dict = read_region_from_bed_file(
        close_chromatin_v5_bed_file)

    gwas_asso_df_h5 = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/test/ai.gemma.gwas_df.h5"
    gwas_asso_df = pd.read_hdf(gwas_asso_df_h5, key="df")

    test_tag = 'CHROMATIN_STATE'
    foreground_tag = 'OPEN_CHROMATIN'
    background_tag = 'CLOSE_CHROMATIN'
    permutation_num = 1000
    test_value = 'LOG_P_VALUE'

    print("Add region tag for each site")
    gwas_asso_df = add_site_tag_by_region(
        gwas_asso_df, test_tag, foreground_tag, open_chromatin_dict, background_tag, close_chromatin_dict)
    gwas_asso_df['LOG_P_VALUE'] = -np.log10(gwas_asso_df['P_VALUE'])

    print("All sites")
    all_maf_observe_ratio, all_maf_perm_ratio_list, all_maf_p_value = site_region_enrichment_test(
        gwas_asso_df, test_tag, foreground_tag, background_tag, test_value, permutation_num, threads=10)

    # maf 0 - 0.5, step 0.1
    maf_result_dict = {}
    for maf in np.arange(0.1, 0.6, 0.1):
        print("MAF: from %.1f to %.1f" % (maf - 0.1, maf))
        gwas_asso_maf_df = gwas_asso_df[(gwas_asso_df['MAF'] >= (
            maf - 0.1)) & (gwas_asso_df['MAF'] < maf)]
        maf_observe_ratio, maf_perm_ratio_list, maf_p_value = site_region_enrichment_test(
            gwas_asso_maf_df, test_tag, foreground_tag, background_tag, test_value, permutation_num, threads=10)
        maf_result_dict[(maf - 0.1, maf)] = (maf_observe_ratio,
                                             maf_perm_ratio_list, maf_p_value)

    for maf in maf_result_dict:
        print(maf, maf_result_dict[maf][0], maf_result_dict[maf][2])
