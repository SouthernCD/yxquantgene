import numpy as np

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
            f.write(sample_name + ' ' + ' '.join(["%.6f" % x for x in dis_matrix[i]]) + '\n')

if __name__ == '__main__':
    matrix_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.ibs_dist.matrix"
    matrix = read_matrix_from_file(matrix_file)

    phylip_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.ibs_dist.phylip"
    sample_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/landraces.id.list"

    from yxutil import read_list_file
    sample_list = read_list_file(sample_file)

    write_matrix_to_phylp_file(matrix, sample_list, phylip_file)