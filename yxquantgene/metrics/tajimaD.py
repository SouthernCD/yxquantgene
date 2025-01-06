import allel
from cyvcf2 import VCF
import numpy as np
import pandas as pd
from yxquantgene.utils.utils import get_chr_len_df
from yxmath.split import split_sequence_to_bins
from yxutil import log_print
import gc


def get_tajima_d_for_genome(vcf_file, reference_genome, window_size, output_h5_file, min_sites=3):
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

            counts, tajima_d_values, windows = get_tajima_d_values(
                genotypes, positions, window_size, start, end, min_sites=min_sites)

            result_df = pd.DataFrame({
                "chr_id": chr_id,
                "start": [i[0] for i in windows],
                "end": [i[1] for i in windows],
                "tajima_d": tajima_d_values,
                "sites": counts
            })

            result_df.to_hdf(output_h5_file, key=chr_id, mode='a', append=True)

            del genotypes
            del positions
            del counts
            del tajima_d_values
            del windows
            del result_df
            gc.collect()

    vcf.close()

    return output_h5_file


def get_tajima_d_for_region(vcf_file, chr_id, start, end, window_size, min_sites=3):
    genotypes, positions = load_genotypes(vcf_file, chr_id, start, end)
    counts, tajima_d_values, windows = get_tajima_d_values(
        genotypes, positions, window_size, start, end, min_sites)

    result_df = pd.DataFrame({
        "chr_id": chr_id,
        "bin": pd.IntervalIndex.from_arrays([i[0] for i in windows], [i[1] for i in windows], closed='both'),
        "tajima_d": tajima_d_values,
        "sites": counts
    })

    return result_df


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


def get_tajima_d_values(genotypes, positions, window_size, start, end, min_sites=3):
    # 将基因型转换为 scikit-allel 支持的格式
    g = allel.GenotypeArray(genotypes)

    # 将基因型数据转换为等位基因数目格式
    ac = g.count_alleles()

    # 计算 Tajima's D
    tajima_d_values, windows, counts = allel.windowed_tajima_d(
        positions, ac, size=window_size, start=start, stop=end, min_sites=min_sites, step=window_size//2)

    # 打印结果
    # for tajima_d, window in zip(tajima_d_values, windows):
    #     print(f"Window {window[0]}-{window[1]}: Tajima's D = {tajima_d}")

    return counts, tajima_d_values, windows


if __name__ == "__main__":
    vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Data/reseq/Sorghum_d8.noduplicates.allChr.snp._markernamesadded_imputed_snpeff.vcf.gz"
    reference_genome = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    chr_id = "Chr01"
    start = 6000000
    end = 7000000
    window_size = 10000

    result_df = get_tajima_d_for_region(
        vcf_file, chr_id, start, end, window_size)

    window_size = 10000
    output_h5_file = f"/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/0.reference/Data/reseq/Sorghum_d8.noduplicates.allChr.snp._markernamesadded_imputed_snpeff.tajimaD.{window_size}.h5"
    get_tajima_d_for_genome(vcf_file, reference_genome,
                            window_size, output_h5_file)

    chr_id = 'Chr01'
    tajimad_df = pd.read_hdf(output_h5_file, key=chr_id)
