import numpy as np
from yxmath.split import split_sequence_to_bins, bin_index_to_range
from yxutil import multiprocess_running
from cyvcf2 import VCF
from yxseq import read_fasta_by_faidx
from itertools import combinations

"""
This module provides functions to calculate Linkage Disequilibrium (LD) matrix.
"""


def calculate_LD(genotype_matrix):
    """
    Calculate LD matrix for a genotype matrix.
    genotype_matrix.shape = (n_variants, n_samples)
    genotype_matrix = np.array([[0, 0, 1, 1],
                                [0, 1, 0, 1],
                                [0, 0, 0, 1],
                                [0, 1, 1, 1],
                                [0, 0, 1, 1]])
    """
    correlation_matrix = np.corrcoef(genotype_matrix)
    LD_matrix = correlation_matrix ** 2
    return LD_matrix


def calculate_LD_by_window(genotype_matrix, window_size=10000, n_jobs=8):
    """
    Calculate LD matrix for a genotype matrix by window, matrix should from the same chromosome.
    genotype_matrix.shape = (n_variants, n_samples)
    genotype_matrix = np.array([[0, 0, 1, 1],
                                [0, 1, 0, 1],
                                [0, 0, 0, 1],
                                [0, 1, 1, 1],
                                [0, 0, 1, 1],
                                [0, 0, 1, 1],
                                [0, 1, 0, 1],
                                [0, 0, 0, 1],
                                [0, 1, 1, 1],
                                [0, 0, 1, 1],
                                [0, 0, 1, 1]])
    """

    n_variants, n_samples = genotype_matrix.shape

    args_dict = {}

    for i, start_tmp, end_tmp in split_sequence_to_bins(n_variants-1, window_size, start=0):
        genotype_matrix_window = genotype_matrix[start_tmp:end_tmp+1, :]

        args_dict[i] = (genotype_matrix_window,)

    mlt_dict = multiprocess_running(calculate_LD, args_dict, n_jobs)

    LD_matrix_dict = {i: mlt_dict[i]['output'] for i in mlt_dict}

    return LD_matrix_dict


def get_chromosome_info(fasta_file):
    fa_dict = read_fasta_by_faidx(fasta_file)
    chrom_len = {}
    for i in fa_dict:
        chrom_len[i] = fa_dict[i].len()
    return chrom_len


def calculate_LD_for_vcf_file(vcf_file, ref_genome_fa, output_prefix, window_size=10000, min_maf=0.1, max_missf=0.5, r2_threshold=0.5):
    vcf = VCF(vcf_file)
    chrom_len_dict = get_chromosome_info(ref_genome_fa)
    chrom_len_dict = {i: chrom_len_dict[i] for i in vcf.seqnames}

    ld_output_file = output_prefix + '.win%d.maf%.2f.miss%.2f.ld' % (
        window_size, min_maf, max_missf)
    var_info_output_file = output_prefix + '.win%d.maf%.2f.miss%.2f.var_info' % (
        window_size, min_maf, max_missf)
    ld_representative_file = output_prefix + '.win%d.maf%.2f.miss%.2f.rq%.2f.ld.nr.rep' % (
        window_size, min_maf, max_missf, r2_threshold)
    ld_remove_file = output_prefix + '.win%d.maf%.2f.miss%.2f.rq%.2f.ld.nr.remove' % (
        window_size, min_maf, max_missf, r2_threshold)

    ld_output_file_handle = open(ld_output_file, 'w')
    var_info_output_file_handle = open(var_info_output_file, 'w')
    ld_representative_file_handle = open(ld_representative_file, 'w')
    ld_remove_file_handle = open(ld_remove_file, 'w')

    num = 0
    for chr_id in chrom_len_dict:
        chrom_len = chrom_len_dict[chr_id]
        for _, s, e in split_sequence_to_bins(chrom_len, window_size, start=1):
            print(f'Processing {chr_id}:{s}-{e}')

            genotypes = []
            var_info_dict = {}
            num = 0
            for v in vcf(f'{chr_id}:{s}-{e}'):
                variant_genotypes = [
                    gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in v.genotypes]
                missf = variant_genotypes.count(-1) / len(variant_genotypes)
                no_miss_sample = len(variant_genotypes) - \
                    variant_genotypes.count(-1)
                maf = min(variant_genotypes.count(0) / no_miss_sample,
                          (variant_genotypes.count(1) + variant_genotypes.count(2)) / no_miss_sample)
                # print((v.ID, v.CHROM, v.POS, v.REF, v.ALT[0], missf, maf))
                if missf <= max_missf and maf >= min_maf:
                    genotypes.append(variant_genotypes)
                    var_info_dict[num] = (
                        v.ID, v.CHROM, v.POS, v.REF, v.ALT[0], missf, maf)
                    num += 1
            genotype_matrix = np.array(genotypes)
            if len(genotype_matrix) == 0:
                continue
            elif len(genotype_matrix) == 1:
                v_ID, v_CHROM, v_POS, v_REF, v_ALT, v_missf, v_maf = var_info_dict[0]
                var_info_output_file_handle.write(
                    f'{v_ID}\t{v_CHROM}\t{v_POS}\t{v_REF}\t{v_ALT}\t{v_missf}\t{v_maf}\n')
                ld_representative_file_handle.write(f'{v_ID}\n')

            elif len(genotype_matrix) > 1:
                LD_matrix = calculate_LD(genotype_matrix)

                linked_pair_list = []
                for i, j in combinations(range(len(genotype_matrix)), 2):
                    v1_ID, v1_CHROM, v1_POS, v1_REF, v1_ALT, v1_missf, v1_maf = var_info_dict[i]
                    v2_ID, v2_CHROM, v2_POS, v2_REF, v2_ALT, v2_missf, v2_maf = var_info_dict[j]
                    ld_output_file_handle.write(
                        f'{v1_ID}\t{v2_ID}\t{LD_matrix[i,j]}\n')
                    if LD_matrix[i, j] >= r2_threshold:
                        linked_pair_list.append(tuple(sorted([v1_ID, v2_ID])))

                for i in var_info_dict:
                    v_ID, v_CHROM, v_POS, v_REF, v_ALT, v_missf, v_maf = var_info_dict[i]
                    var_info_output_file_handle.write(
                        f'{v_ID}\t{v_CHROM}\t{v_POS}\t{v_REF}\t{v_ALT}\t{v_missf}\t{v_maf}\n')

                representative_variants = []
                remove_variants = []

                for i in sorted(var_info_dict, key=lambda x: var_info_dict[x][6], reverse=True):
                    v_ID, v_CHROM, v_POS, v_REF, v_ALT, v_missf, v_maf = var_info_dict[i]
                    linked_variants = get_linked_variants(
                        v_ID, linked_pair_list)
                    if len(set(representative_variants) & set(linked_variants)) == 0:
                        representative_variants.append(v_ID)
                    else:
                        remove_variants.append(v_ID)

                var_info_dict = {
                    var_info_dict[i][0]: var_info_dict[i] for i in var_info_dict}

                representative_variants = list(set(representative_variants))
                representative_variants = sorted(
                    representative_variants, key=lambda x: var_info_dict[x][2], reverse=False)
                remove_variants = list(set(remove_variants))
                remove_variants = sorted(
                    remove_variants, key=lambda x: var_info_dict[x][2], reverse=False)

                for i in representative_variants:
                    ld_representative_file_handle.write(f'{i}\n')
                for i in remove_variants:
                    ld_remove_file_handle.write(f'{i}\n')

            num += 1

            # if num > 10:
            #     break

    ld_output_file_handle.close()
    var_info_output_file_handle.close()
    ld_representative_file_handle.close()
    ld_remove_file_handle.close()


def get_linked_variants(v_ID, linked_pair_list):
    """
    v_ID = '1'
    linked_pair_list = [('1', '2'), ('3', '1'), ('2', '3'), ('2', '4'), ('3', '4'), ('1', '1')]
    """
    linked_variants = []
    for v1, v2 in linked_pair_list:
        if v1 == v_ID:
            linked_variants.append(v2)
        elif v2 == v_ID:
            linked_variants.append(v1)
    linked_variants = list(set(linked_variants))
    if v_ID in linked_variants:
        linked_variants.remove(v_ID)
    return linked_variants



if __name__ == '__main__':
    ref_genome_fa = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa'
    vcf_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.vcf.gz'

    output_prefix = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp'

    calculate_LD_for_vcf_file(vcf_file, ref_genome_fa, output_prefix,
                              window_size=10000, min_maf=0.1, max_missf=0.5, r2_threshold=0.5)
