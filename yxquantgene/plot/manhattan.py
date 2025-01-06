import matplotlib.pyplot as plt
import numpy as np


def local_general_plotter(data_df, chr_id, start, end, ax=None, style='scatter', c1="#617EB8", ylabel=""):
    """
    support style: scatter, line, shaded

    data_df: DataFrame with columns: CHROM, POS, VALUE
    chr_id: chromosome id
    start: start position
    end: end position
    style: scatter, line, shaded
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    else:
        fig = ax.get_figure()

    data_df = data_df.dropna()
    data_df['CHROM'] = data_df['CHROM'].astype(str)
    data_df['POS'] = data_df['POS'].astype(int)
    data_df['VALUE'] = data_df['VALUE'].astype(float)

    data_df = data_df.loc[(data_df['CHROM'] == chr_id) & (
        data_df['POS'] >= start) & (data_df['POS'] <= end)]

    if data_df.shape[0] == 0:
        return

    if style == 'scatter':
        ax.scatter(data_df['POS'], data_df['VALUE'], s=10, c=c1)
    elif style == 'line':
        ax.plot(data_df['POS'], data_df['VALUE'], c=c1)
    elif style == 'shaded':
        ax.fill_between(data_df['POS'], data_df['VALUE'], color=c1, alpha=0.5)

    ax.set_xlim(start, end)
    ax.set_ylim(np.min(data_df['VALUE']), np.max(data_df['VALUE'])*1.1)

    ax.set_xlabel(f'{chr_id}: {start}-{end}')
    ax.set_ylabel(ylabel)

    if ax is None:
        plt.show()


def local_gwas_plot(gwas_df, chr_id, start, end, ax=None, pval_lim=1e-20):
    """
    绘制 Manhattan 图
    gwas_df should have columns: CHROM, POS, P_VALUE, FDR
    chr_list: list of chromosome names
    chr_length_dict: dict of chromosome length
    """
    # data preparation
    gwas_df.loc[gwas_df['P_VALUE'] <= pval_lim, 'P_VALUE'] = pval_lim
    gwas_df.loc[gwas_df['FDR'] == 0, 'FDR'] = pval_lim

    data_df = gwas_df.copy()
    data_df['VALUE'] = -np.log10(data_df['P_VALUE'])

    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    else:
        fig = ax.get_figure()

    local_general_plotter(data_df, chr_id, start, end, ax=ax,
                          style='scatter', ylabel='$-log_{10}(p)$')

    if ax is None:
        plt.show()


def genomewide_general_plotter(data_df, chr_length_df, ax=None, style='scatter', c1="#617EB8", c2="#84B8D0", min_chr_length=1e6, max_y=None, ylabel=""):
    """
    support style: scatter, line, shaded

    data_df: DataFrame with columns: CHROM, POS, VALUE
    chr_length_df: DataFrame with columns: chr, length. Chrom rank in chr_length_df will be used to determine the plot order
    style: scatter, line, shaded
    """

    chr_length_df = chr_length_df.loc[chr_length_df['LEN'] > min_chr_length]

    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    else:
        fig = ax.get_figure()

    chr_coord_dict = {}
    total_length = 0
    chr_list = chr_length_df['CHROM'].tolist()
    for i, r in chr_length_df.iterrows():
        chr_coord_dict[r['CHROM']] = total_length
        total_length += r['LEN']

    # 设置ax
    data_df = data_df.dropna()
    data_df['CHROM'] = data_df['CHROM'].astype(str)
    data_df['POS'] = data_df['POS'].astype(int)
    data_df['VALUE'] = data_df['VALUE'].astype(float)

    ax.set_xlim(0 - total_length*0.05, total_length*1.05)
    max_y = np.max(data_df['VALUE']) if max_y is None else max_y
    ax.set_ylim(0, max_y*1.1)

    # 设置 x 轴和 y 轴的标签
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylabel)

    # 隐藏坐标轴的上边框和右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 绘图
    c1 = "#617EB8"
    c2 = "#84B8D0"
    n = 0
    chr_ticks = []  # 用于存储每个染色体中点的位置
    chr_labels = []  # 用于存储染色体的编号
    for chr_id in chr_list:
        c = c1 if n % 2 == 0 else c2
        chr_df = data_df.loc[data_df['CHROM'] == chr_id]
        if style == 'scatter':
            ax.scatter((chr_df['POS'] + chr_df['CHROM'].map(chr_coord_dict)
                        ).values, chr_df['VALUE'].values, s=10, c=c)
        elif style == 'line':
            ax.plot((chr_df['POS'] + chr_df['CHROM'].map(chr_coord_dict)
                     ).values, chr_df['VALUE'].values, c=c)
        elif style == 'shaded':
            ax.fill_between((chr_df['POS'] + chr_df['CHROM'].map(chr_coord_dict)
                             ).values, chr_df['VALUE'].values, color=c, alpha=0.5)
        # 计算染色体的中点位置并添加到列表中
        chr_ticks.append(
            (chr_df['POS'] + chr_df['CHROM'].map(chr_coord_dict)).mean())
        chr_labels.append(chr_id)
        n += 1

    # 设置 x 轴的刻度标签
    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(chr_labels)

    if ax is None:
        plt.show()


def genomewide_gwas_plot(gwas_df, chr_length_df, threshold_fdr=0.05, ax=None, pval_lim=1e-20):
    """
    绘制 Manhattan 图
    gwas_df should have columns: CHROM, POS, P_VALUE, FDR
    chr_list: list of chromosome names
    chr_length_dict: dict of chromosome length
    """
    # data preparation
    gwas_df.loc[gwas_df['P_VALUE'] <= pval_lim, 'P_VALUE'] = pval_lim
    gwas_df.loc[gwas_df['FDR'] == 0, 'FDR'] = pval_lim

    data_df = gwas_df.copy()
    data_df['VALUE'] = -np.log10(data_df['P_VALUE'])

    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    else:
        fig = ax.get_figure()

    genomewide_general_plotter(data_df, chr_length_df, ax=ax,
                               style='scatter', ylabel='$-log_{10}(p)$')

    # 找到qval阈值对应的-log10(pval)
    threshold = - \
        np.log10(gwas_df.loc[gwas_df['FDR'] < threshold_fdr, 'P_VALUE'].max())
    ax.axhline(y=threshold, color='r', linestyle='--',
               label='FDR == %.2f' % threshold_fdr)

    if ax is None:
        plt.show()


if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    from yxquantgene.utils.utils import get_chr_len_df

    # whole genome association study plot

    genome_fasta_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa"
    gwas_df_h5 = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/all_landraces_gemma/ai/snp/gemma_output.assoc.h5"

    gwas_df = pd.read_hdf(gwas_df_h5, 'df')
    gwas_df = gwas_df.rename(
        columns={'CHROM': 'CHROM', 'POS': 'POS', 'pvalue': 'P_VALUE', 'fdr_bh': 'FDR'})
    chr_length_df = get_chr_len_df(genome_fasta_file)

    genomewide_gwas_plot(gwas_df, chr_length_df,
                         threshold_fdr=0.05, ax=None, pval_lim=1e-20)

    # local association study plot
