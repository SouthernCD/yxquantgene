import numpy as np
from yxseq import read_fasta_by_faidx, read_gff_file
from yxutil import have_file
import pandas as pd


def get_chr_len_df(fasta_file):
    chrom_len_csv = f'{fasta_file}.ChrLen.csv'
    if have_file(chrom_len_csv):
        chrom_len_df = pd.read_csv(chrom_len_csv)
    else:
        fa_dict = read_fasta_by_faidx(fasta_file)
        chrom_len = {}
        for i in fa_dict:
            chrom_len[i] = fa_dict[i].len()
        chrom_len_df = pd.DataFrame(chrom_len.items(), columns=['CHROM', 'LEN'])
        chrom_len_df.to_csv(chrom_len_csv, index=False)
    return chrom_len_df


def get_gene_range_df(gff_file):
    gene_range_csv = f'{gff_file}.GeneRange.csv'
    if have_file(gene_range_csv):
        gene_range_df = pd.read_csv(gene_range_csv)
    else:
        gene_dict = read_gff_file(gff_file)['gene']
        gene_id_list = list(gene_dict.keys())
        gene_range_df = pd.DataFrame()
        gene_range_df['GENE'] = gene_id_list
        gene_range_df['CHROM'] = [gene_dict[i].chr_id for i in gene_id_list]
        gene_range_df['START'] = [gene_dict[i].start for i in gene_id_list]
        gene_range_df['END'] = [gene_dict[i].end for i in gene_id_list]
        gene_range_df['STRAND'] = [gene_dict[i].strand for i in gene_id_list]
        gene_range_df.to_csv(gene_range_csv, index=False)
    return gene_range_df


def write_matrix_to_file(matrix, file_name):
    """
    将矩阵写入文件。
    """
    np.savetxt(file_name, matrix, fmt='%.6f', delimiter=', ')


def read_matrix_from_file(file_name):
    """
    从文件中读取矩阵。
    """
    matrix = np.loadtxt(file_name, delimiter=',')
    return matrix


def write_matrix_to_phylp_file(dis_matrix, sample_names, file_name):
    """
    将距离矩阵写入PHYLIP格式文件。
    """
    with open(file_name, 'w') as f:
        f.write(str(len(sample_names)) + '\n')
        for i, sample_name in enumerate(sample_names):
            f.write(sample_name + ' ' +
                    ' '.join(["%.6f" % x for x in dis_matrix[i]]) + '\n')


if __name__ == '__main__':
    matrix_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.ibs_dist.matrix"
    matrix = read_matrix_from_file(matrix_file)

    phylip_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.ibs_dist.phylip"
    sample_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/landraces.id.list"

    from yxutil import read_list_file
    sample_list = read_list_file(sample_file)

    write_matrix_to_phylp_file(matrix, sample_list, phylip_file)
