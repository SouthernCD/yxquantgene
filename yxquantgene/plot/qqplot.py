import matplotlib.pyplot as plt
import numpy as np

def quantile_quantile_plot(p_value, title=None, ax=None, filter_zero_point=True):
    p_value = np.array(p_value)

    if ax is None:
        fig, ax = plt.subplots()

    if filter_zero_point:
        # 计算 -np.log10(p_value)
        log_p_value = -np.log10(p_value)
        # 找出 -np.log10(p_value) 不等于 0 的值
        mask = log_p_value != 0
        # 使用 mask 来筛选 p_value
        p_value = p_value[mask]

    # 计算理论分位数
    n = len(p_value)

    if n == 0:
        ax.plot([0, 1], [0, 1], c="#FB0324")
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
    else:
        expected = -np.log10(np.linspace(1/n, 1, n))

        max_value = max(np.max(expected), np.max(-np.log10(p_value)))

        # 绘图
        ax.scatter(expected, -np.log10(np.sort(p_value)), s=20, c="#000000")

        lim = (0 - max_value*0.1, max_value*1.1)
        ax.plot(lim, lim, c="#FB0324")

        # 设置 x 轴和 y 轴的范围
        ax.set_xlim(0 - np.max(expected)*0.1, np.max(expected)*1.1)
        ax.set_ylim(0 - np.max(-np.log10(p_value))*0.1, np.max(-np.log10(p_value))*1.1)

    # 设置 x 轴和 y 轴的标签
    ax.set_xlabel('Expected $-log_{10}(p)$')
    ax.set_ylabel('Observed $-log_{10}(p)$')

    # 设置标题
    if title is not None:
        ax.set_title(title)

    if ax is None:
        plt.show()