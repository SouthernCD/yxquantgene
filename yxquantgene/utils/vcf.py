from cyvcf2 import VCF, Writer
from yxutil import cmd_run
import numpy as np


def extract_subvcf_by_varIDs(input_vcf_file, varID_list, output_vcf_file):
    """
    Extract a subset of variants from a VCF file by a list of variant IDs.
    """
    varID_list = set(varID_list)

    vcf_reader = VCF(input_vcf_file)
    vcf_writer = Writer(output_vcf_file, vcf_reader)

    for record in vcf_reader:
        if record.ID in varID_list:
            vcf_writer.write_record(record)

    vcf_writer.close()
    vcf_reader.close()

    cmd_run(f'bgzip {output_vcf_file}')
    cmd_run(f'tabix -f -p vcf {output_vcf_file}.gz')

    return output_vcf_file + '.gz'


def get_genotype_matrix_from_vcf(vcf_file, chr_id=None, start=None, end=None):
    """
    Get genotype matrix from a VCF file.
    """

    vcf = VCF(vcf_file)
    genotype_matrix = []

    if chr_id is not None and start is not None and end is not None:
        for record in vcf(f'{chr_id}:{start}-{end}'):
            variant_genotypes = [
                gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
            genotype_matrix.append(variant_genotypes)
    else:
        for record in vcf:
            variant_genotypes = [
                gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
            genotype_matrix.append(variant_genotypes)

    genotype_matrix = np.array(genotype_matrix)

    return genotype_matrix


if __name__ == '__main__':
    from yxutil import read_list_file

    varID_list_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep"
    varID_list = read_list_file(varID_list_file)

    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.vcf.gz"
    output_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.vcf"

    output_vcf_file = extract_subvcf_by_varIDs(
        input_vcf_file, varID_list, output_vcf_file)

    genotype_matrix = get_genotype_matrix_from_vcf(output_vcf_file)
