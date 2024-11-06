import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(manhattan_df, chr_list, chr_length_dict, threshold_qval=0.05 ,ax=None):
    """
    绘制 Manhattan 图
    manhattan_df should have columns: chr, pos, pval, qval
    chr_list: list of chromosome names
    chr_length_dict: dict of chromosome length
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(16,9))
    else:
        fig = ax.get_figure()

    chr_coord_dict = {}
    total_length = 0
    for i in chr_list:
        chr_coord_dict[i] = total_length
        total_length += chr_length_dict[i]

    # 设置ax
    ax.set_xlim(0 - total_length*0.05, total_length*1.05)
    max_y = np.max(-np.log10(manhattan_df.pval))
    ax.set_ylim(0, max_y*1.1)

    # 设置 x 轴和 y 轴的标签
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('$-log_{10}(p)$')

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
        chr_df = manhattan_df.loc[manhattan_df['chr'] == chr_id]
        ax.scatter(chr_df.pos + chr_df.chr.map(chr_coord_dict), -np.log10(chr_df.pval), s=10, c=c)
        # 计算染色体的中点位置并添加到列表中
        chr_ticks.append((chr_df.pos + chr_df.chr.map(chr_coord_dict)).mean())
        chr_labels.append(chr_id)
        n += 1

    # 找到qval阈值对应的-log10(pval)
    threshold = -np.log10(manhattan_df.loc[manhattan_df.qval < threshold_qval, 'pval'].max())
    ax.axhline(y=threshold, color='r', linestyle='--', label='FDR == %.2f' % threshold_qval)

    # 设置 x 轴的刻度标签
    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(chr_labels)

    if ax is None:
        plt.show()