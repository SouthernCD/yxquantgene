from sklearn.decomposition import PCA
from yxquantgene.utils.vcf import get_genotype_matrix_from_vcf, get_sample_list_from_vcf
import pandas as pd
from yxutil import log_print


def get_PCA_matrix(input_matrix, n_components=2):
    """
    input_matrix: 输入矩阵，行为变异，列为样品
    """
    pca = PCA(n_components=n_components)
    pca_matrix = pca.fit_transform(input_matrix.T) # 注意这里的转置，
    explained_variance_ratio = pca.explained_variance_ratio_  # 解释方差比例

    return pca_matrix, explained_variance_ratio


def PCA_analysis(input_matrix_df, n_components=2):
    """
    input_matrix_df: 输入dataframe, 行为变异，列为样品
    """

    pca_matrix, explained_variance_ratio = get_PCA_matrix(
        input_matrix_df.values, n_components)

    pca_df = pd.DataFrame(pca_matrix, columns=[
                          f'PC{i+1}' for i in range(n_components)], index=input_matrix_df.index)
    explained_variance_ratio_df = pd.DataFrame(explained_variance_ratio, index=[
                                               f'PC{i+1}' for i in range(n_components)], columns=['Explained Variance Ratio'])

    return pca_df, explained_variance_ratio_df


def PCA_for_vcf_file(input_vcf_file, n_components=2):
    """
    从vcf文件中读取基因型矩阵（行为变异，列为样品），计算PCA矩阵
    """
    log_print('Reading genotype matrix from VCF file...')
    genotype_matrix = get_genotype_matrix_from_vcf(input_vcf_file)
    log_print('Calculating PCA matrix...')
    pca_matrix, explained_variance_ratio = get_PCA_matrix(
        genotype_matrix, n_components)
    log_print('PCA matrix calculated.')
    pca_df = pd.DataFrame(pca_matrix, columns=[
        f'PC{i+1}' for i in range(n_components)], index=get_sample_list_from_vcf(input_vcf_file))
    explained_variance_ratio_df = pd.DataFrame(explained_variance_ratio, index=[
        f'PC{i+1}' for i in range(n_components)], columns=['Explained Variance Ratio'])

    return pca_df, explained_variance_ratio_df

if __name__ == '__main__':
    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/ld.win10000.maf0.01.miss0.50.rq0.50.ld.nr.rep.vcf.gz"
    pca_df, explained_variance_ratio_df = PCA_for_vcf_file(input_vcf_file)
    print(pca_df)
    print(explained_variance_ratio_df)
    pca_df.to_csv('pca.csv')
    explained_variance_ratio_df.to_csv('explained_variance_ratio.csv')