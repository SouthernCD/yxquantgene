import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

def allele_corr_plot(input_df, ref_allele=None, alt_allele=None, phased=False, trendline=True, y_label='Value', title='Correlation', save_file=None, ax=None):
    """
    input_df: pandas.DataFrame
    should have columns: genotype, and value columns
    genotype: genotype of the site (like '0|0', '0|1', '1|1'), 0 for reference allele, 1 for alternative allele
    value: value to be plotted as y-axis

    statistic = 'pearson', 'spearman' or 'kruskal'
    """
    input_df = input_df[input_df['value'].notnull()]
    input_df = input_df[input_df['genotype'].notnull()]

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 6))

    if phased:
        x_label_list = ['0|0', '0|1', '1|0', '1|1']
    else:
        x_label_list = ['0/0', '0/1', '1/1']
        ref_df = input_df.loc[input_df['genotype'].isin(['0|0', '0/0'])].copy()
        ref_df.loc[:, 'genotype'] = '0/0'
        hap_df = input_df.loc[input_df['genotype'].isin(
            ['0|1', '1|0', '1/0', '0/1'])].copy()
        hap_df.loc[:, 'genotype'] = '0/1'
        alt_df = input_df.loc[input_df['genotype'].isin(['1|1', '1/1'])].copy()
        alt_df.loc[:, 'genotype'] = '1/1'
        input_df = pd.concat([ref_df, hap_df, alt_df], ignore_index=True)

    sns.boxplot(x='genotype', y='value', data=input_df, showfliers=False,
                order=x_label_list, width=0.3, color='#BEDDFD', ax=ax)
    sns.stripplot(x='genotype', y='value', data=input_df, color='#258FF8',
                  jitter=0.2, order=x_label_list, size=3, alpha=0.5, ax=ax)

    # x = 0
    # x_label_list = []
    # x_list = []
    # y_list = []

    # df = input_df[input_df['genotype'].isin(['0|0', '0/0'])]
    # ax.scatter(np.full(len(df['value']), x), df['value'],
    #             color='#3B75AF', label='0|0' if phased else '0/0')
    # x_label_list.append('0|0' if phased else '0/0')
    # x_list.extend(list(np.full(len(df['value']), x)))
    # y_list.extend(list(df['value']))
    # x += 1

    # if phased:
    #     df = input_df[input_df['genotype'].isin(['0|1'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='0|1')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('0|1')
    #     x += 1

    #     df = input_df[input_df['genotype'].isin(['1|0'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='1|0')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('1|0')
    #     x += 1

    # else:
    #     df = input_df[input_df['genotype'].isin(['0|1', '1|0', '1/0', '0/1'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='0/1')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('0/1')
    #     x += 1

    # df = input_df[input_df['genotype'].isin(['1|1', '1/1'])]
    # ax.scatter(np.full(len(df['value']), x), df['value'],
    #             color='#3B75AF', label='1|1' if phased else '1/1')
    # x_list.extend(list(np.full(len(df['value']), x)))
    # y_list.extend(list(df['value']))
    # x_label_list.append('1|1' if phased else '1/1')
    # x += 1

    # ax.set_xticks(range(x))
    # ax.set_xticklabels(x_label_list)

    ax.set_xlim(-0.5, len(x_label_list)-0.5)

    if ref_allele is not None and alt_allele is not None:
        x_label = "Genotype (Ref: {0}, Alt: {1})".format(
            ref_allele, alt_allele)
    else:
        x_label = "Genotype"

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if not phased:
        # 拟合数据并绘制趋势线
        if trendline:
            x_list = input_df['genotype']
            x_list.replace('0/0', 0, inplace=True)
            x_list.replace('0/1', 1, inplace=True)
            x_list.replace('1/1', 2, inplace=True)
            y_list = input_df['value']

            try:
                z = np.polyfit(x_list, y_list, 1)
                p = np.poly1d(z)
                ax.plot(x_list, p(x_list), "-", c='#E87D85')
            except:
                pass

        p = pearsonr(x_list, y_list)
        subtitile = "Pearson correlation: %.4e, p-value: %.4e" % (
            p[0], p[1])

        ax.text(0.5, 1.12, title, transform=ax.transAxes,
                ha="center", va="center", fontsize=18)
        ax.text(0.5, 1.05, subtitile, transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
    else:
        ax.text(0.5, 1.05, title, transform=ax.transAxes,
                ha="center", va="center", fontsize=18)

    if ax is None:
        plt.show()

    if save_file is not None:
        fig.savefig(save_file, format='pdf', facecolor='none',
                    edgecolor='none', bbox_inches='tight')

    return ax, p[0], p[1]