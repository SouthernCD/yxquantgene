import allel
from cyvcf2 import VCF
import numpy as np
import pandas as pd
from yxquantgene.utils.utils import get_chr_len_df
from yxmath.split import split_sequence_to_bins
from yxutil import log_print
import gc

# 读取 VCF 文件


def load_genotypes(vcf_file, chr_id, start, end):
    region = f"{chr_id}:{start}-{end}"

    # 初始化 VCF 对象
    vcf = VCF(vcf_file)

    # 使用 fetch 方法读取特定区段
    positions = []
    genotypes = []

    for variant in vcf(region):
        positions.append(variant.POS)
        genotypes.append(variant.genotypes)

    # 转换为 NumPy 格式，scikit-allel 需要该格式
    genotypes = np.array(genotypes)[:, :, :2]  # 去掉相位信息 (只保留两等位基因)
    positions = np.array(positions)

    vcf.close()

    return genotypes, positions


def windowed_r_squared(pos, gn, size=None, start=None, stop=None, step=None,
                       windows=None, fill=np.nan, percentile=50):
    """Summarise linkage disequilibrium in windows over a single
    chromosome/contig.

    Parameters
    ----------
    pos : array_like, int, shape (n_items,)
        The item positions in ascending order, using 1-based coordinates..
    gn : array_like, int8, shape (n_variants, n_samples)
        Diploid genotypes at biallelic variants, coded as the number of
        alternate alleles per call (i.e., 0 = hom ref, 1 = het, 2 = hom alt).
    size : int, optional
        The window size (number of bases).
    start : int, optional
        The position at which to start (1-based).
    stop : int, optional
        The position at which to stop (1-based).
    step : int, optional
        The distance between start positions of windows. If not given,
        defaults to the window size, i.e., non-overlapping windows.
    windows : array_like, int, shape (n_windows, 2), optional
        Manually specify the windows to use as a sequence of (window_start,
        window_stop) positions, using 1-based coordinates. Overrides the
        size/start/stop/step parameters.
    fill : object, optional
        The value to use where a window is empty, i.e., contains no items.
    percentile : int or sequence of ints, optional
        The percentile or percentiles to calculate within each window.

    Returns
    -------
    out : ndarray, shape (n_windows,)
        The value of the statistic for each window.
    windows : ndarray, int, shape (n_windows, 2)
        The windows used, as an array of (window_start, window_stop) positions,
        using 1-based coordinates.
    counts : ndarray, int, shape (n_windows,)
        The number of items in each window.

    Notes
    -----
    Linkage disequilibrium (r**2) is calculated using the method of Rogers
    and Huff (2008).

    See Also
    --------

    allel.stats.window.windowed_statistic

    """

    # define the statistic function
    if isinstance(percentile, (list, tuple)):
        fill = [fill for _ in percentile]

        def statistic(gnw):
            r_squared = allel.rogers_huff_r(gnw) ** 2
            return [np.percentile(r_squared, p) for p in percentile]

    else:
        def statistic(gnw):
            r_squared = allel.rogers_huff_r(gnw) ** 2
            # print(r_squared)
            if isinstance(r_squared, np.float64):
                return r_squared
            if isinstance(r_squared, np.ndarray) and len(r_squared) == 0:
                return fill
            return np.percentile(r_squared, percentile)

    return allel.windowed_statistic(pos, gn, statistic, size, start=start,
                                    stop=stop, step=step, windows=windows, fill=fill)


def get_win_stat_for_genome(vcf_file, reference_genome, window_size, output_h5_file, min_sites=3):
    chr_len_df = get_chr_len_df(reference_genome)
    split_lenght = 5000000 if window_size < 5000000 else (window_size * 10)

    vcf = VCF(vcf_file)

    for chr_id in chr_len_df['chr_id'].to_list():
        chr_length = chr_len_df[chr_len_df['chr_id']
                                == chr_id]['len'].values[0]

        for i, start, end in split_sequence_to_bins(chr_length, split_lenght, start=1):
            region = f"{chr_id}:{start}-{end}"
            log_print(f"Processing {region}")

            # 使用 fetch 方法读取特定区段
            positions = []
            genotypes = []

            for variant in vcf(region):
                positions.append(variant.POS)
                genotypes.append(variant.genotypes)

            if len(positions) == 0:
                continue

            # 转换为 NumPy 格式，scikit-allel 需要该格式
            genotypes = np.array(genotypes)[:, :, :2]  # 去掉相位信息 (只保留两等位基因)
            positions = np.array(positions)

            # 进行计算
            # 将基因型转换为 scikit-allel 支持的格式
            g = allel.GenotypeArray(genotypes)
            ac = g.count_alleles()
            af = ac.to_frequencies()

            # 计算 Tajima's D
            tajima_d_values, windows, counts = allel.windowed_tajima_d(
                positions, ac, size=window_size, start=start, stop=end, min_sites=min_sites, step=window_size//2)

            # 计算核苷酸多样性 (π)
            pi_values, windows, n_bases, counts = allel.windowed_diversity(
                positions, ac, size=window_size, start=start, stop=end, step=window_size//2)

            # 计算窗口 Watterson's theta
            windowed_watterson_theta, windows, n_bases, counts = allel.windowed_watterson_theta(
                positions, ac, size=window_size, start=start, stop=end, step=window_size//2)

            # 计算杂合度
            ho = allel.heterozygosity_observed(g)
            win_ho, windows, counts = allel.windowed_statistic(
                positions, ho, statistic=np.mean, size=window_size, start=start, stop=end, step=window_size//2)
            he = allel.heterozygosity_expected(af, ploidy=g.shape[-1])
            win_he, windows, counts = allel.windowed_statistic(
                positions, he, statistic=np.mean, size=window_size, start=start, stop=end, step=window_size//2)

            # 计算LD
            gn = g.to_n_alt()
            r2_values, windows, counts = windowed_r_squared(
                positions, gn, size=window_size, start=start, stop=end, step=window_size//2)

            # 保存结果
            result_df = pd.DataFrame({
                "chr_id": chr_id,
                "start": [i[0] for i in windows],
                "end": [i[1] for i in windows],
                "sites": counts,
                "tajima_d": tajima_d_values,
                "pi": pi_values,
                "watterson_theta": windowed_watterson_theta,
                "heterozygosity_observed": win_ho,
                "heterozygosity_expected": win_he,
                "r2": r2_values
            })

            result_df.to_hdf(output_h5_file, key=chr_id, mode='a', append=True)

            del genotypes
            del positions
            del g
            del ac
            del result_df
            gc.collect()

    vcf.close()

    return output_h5_file


if __name__ == '__main__':
    vcf_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces_snp.vcf.gz'
    reference_genome = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    window_size = 100
    output_h5_file = vcf_file + ".win_stat.%d.h5" % window_size
    get_win_stat_for_genome(vcf_file, reference_genome,
                            window_size, output_h5_file, min_sites=3)
    window_size = 1000
    output_h5_file = vcf_file + ".win_stat.%d.h5" % window_size
    get_win_stat_for_genome(vcf_file, reference_genome,
                            window_size, output_h5_file, min_sites=3)
    window_size = 10000
    output_h5_file = vcf_file + ".win_stat.%d.h5" % window_size
    get_win_stat_for_genome(vcf_file, reference_genome,
                            window_size, output_h5_file, min_sites=3)
