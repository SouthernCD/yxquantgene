from sklearn.cluster import KMeans
import pandas as pd
from yxquantgene.utils.vcf import get_genotype_matrix_from_vcf, get_sample_list_from_vcf

def cluster_samples_from_vcf(vcf_file, n_clusters=3, random_state=42):
    """
    从VCF文件中读取基因型数据，并对样品进行KMeans聚类。

    参数:
    - vcf_file (str): VCF文件路径。
    - n_clusters (int): 聚类的簇数。
    - random_state (int): 随机种子，确保结果可重复。

    返回:
    - clustered_samples (DataFrame): 包含样品名和簇标签的DataFrame。
    """
    # 读取基因型矩阵和样品名称
    genotype_matrix = get_genotype_matrix_from_vcf(vcf_file)
    sample_list = get_sample_list_from_vcf(vcf_file)

    # 使用现有的 kmeans_clustering 函数进行聚类
    _, labels, clustered_data = kmeans_clustering(
        data=genotype_matrix, n_clusters=n_clusters, random_state=random_state
    )

    # 添加样品名称
    clustered_data['Sample'] = sample_list

    # 调整列顺序
    clustered_data = clustered_data[['Sample', 'Cluster']]

    return clustered_data


def kmeans_clustering(data, n_clusters=3, random_state=42):
    """
    对数据进行KMeans聚类。
    
    参数:
    - data (array-like 或 DataFrame): 数据集，形状为 (n_samples, n_features)。
    - n_clusters (int): 聚类的簇数。
    - random_state (int): 随机种子，确保结果可重复。
    
    返回:
    - cluster_centers (ndarray): 聚类中心，形状为 (n_clusters, n_features)。
    - labels (ndarray): 每个样本的簇标签，形状为 (n_samples,)。
    - clustered_data (DataFrame): 包含原始数据和簇标签的DataFrame。
    """
    # 初始化KMeans模型
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
    
    # 拟合模型并预测标签
    labels = kmeans.fit_predict(data)
    
    # 获取聚类中心
    cluster_centers = kmeans.cluster_centers_
    
    # 如果输入是DataFrame，将结果附加到DataFrame中
    if isinstance(data, pd.DataFrame):
        clustered_data = data.copy()
        clustered_data['Cluster'] = labels
    else:
        clustered_data = pd.DataFrame(data, columns=[f'Feature_{i+1}' for i in range(data.shape[1])])
        clustered_data['Cluster'] = labels
    
    return cluster_centers, labels, clustered_data

# 示例用法
if __name__ == "__main__":
    pass