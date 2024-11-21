from cyvcf2 import VCF, Writer
from yxquantgene.utils.utils import get_chromosome_info
from yxutil import cmd_run, log_print, multiprocess_running
import gc
import numpy as np
import pandas as pd

def extract_subvcf_by_samples(input_vcf_file, sample_list, output_vcf_file):
    tmp_sample_list_file = 'sample.id.list'
    with open(tmp_sample_list_file, 'w') as f:
        for sample in sample_list:
            f.write(sample + '\n')

    if not output_vcf_file.endswith('.gz'):
        output_vcf_file = output_vcf_file + '.gz'

    cmd_run(f"bcftools view --threads 20 -c 1 -O z -o {output_vcf_file} -S sample.id.list {input_vcf_file}")
    cmd_run(f'tabix -f -p vcf {output_vcf_file}')
    cmd_run(f'rm {tmp_sample_list_file}')


def extract_subvcf_by_varIDs(input_vcf_file, varID_list, output_vcf_file):
    """
    Extract a subset of variants from a VCF file by a list of variant IDs.
    """
    if output_vcf_file.endswith('.gz'):
        output_vcf_file = output_vcf_file[:-3]

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
    elif chr_id is not None and start is None and end is None:
        for record in vcf(f'{chr_id}'):
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


def get_genotype_df_from_vcf(vcf_file, chr_id=None, start=None, end=None):
    """
    Get genotype df from a VCF file.
    """

    vcf = VCF(vcf_file)
    sample_list = vcf.samples
    genotype_matrix = []
    varID_list = []

    num = 0
    if chr_id is not None and start is not None and end is not None:
        for record in vcf(f'{chr_id}:{start}-{end}'):
            variant_genotypes = [
                gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
            genotype_matrix.append(variant_genotypes)
            varID_list.append(record.ID)
            num += 1
            if num % 5000 == 0:
                print(f'Processed {num} variants.')
    else:
        for record in vcf:
            variant_genotypes = [
                gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
            genotype_matrix.append(variant_genotypes)
            varID_list.append(record.ID)
            num += 1
            if num % 5000 == 0:
                print(f'Processed {num} variants.')

    vcf.close()

    if len(varID_list) == 0:
        return None
    else:
        genotype_matrix = np.array(genotype_matrix)
        genotype_df = pd.DataFrame(
            genotype_matrix, columns=sample_list, index=varID_list)

        return genotype_df


def get_sample_list_from_vcf(vcf_file):
    """
    Get sample list from a VCF file.
    """

    vcf = VCF(vcf_file)
    sample_list = vcf.samples

    return sample_list


def get_varID_list_from_vcf(vcf_file, chr_id=None, start=None, end=None):
    """
    Get variant ID list from a VCF file.
    """

    if chr_id is None:
        vcf = VCF(vcf_file)
        varID_list = [record.ID for record in vcf]
        vcf.close()
    else:
        vcf = VCF(vcf_file)
        varID_list = [record.ID for record in vcf(f'{chr_id}:{start}-{end}')]
        vcf.close()

    return varID_list


# variant statistics table
def get_var_stat(genotype_matrix):
    # Count the number of each element in the matrix
    counts = np.apply_along_axis(lambda x: np.histogram(
        x, bins=[-1.5, -0.5, 0.5, 1.5, 2.5])[0], 1, genotype_matrix)

    # Calculate the number of each type of element
    mis_num = counts[:, 0]
    hom_ref_num = counts[:, 1]
    het_num = counts[:, 2]
    hom_alt_num = counts[:, 3]

    maf = (het_num + 2*hom_alt_num)/(2*(hom_ref_num + het_num + hom_alt_num))
    maf = np.minimum(maf, 1-maf)

    het = het_num/(hom_ref_num + het_num + hom_alt_num)
    mis = mis_num/(hom_ref_num + het_num + hom_alt_num + hom_alt_num)

    return mis, maf, het, hom_ref_num, het_num, hom_alt_num


# def get_var_stat_chunk(k, chunk_size, genotype_matrix):
#     k_end = min(k + chunk_size, genotype_matrix.shape[0])
#     genotype_matrix_chunk = genotype_matrix[k:k_end]
#     mis, maf, het, hom_ref_num, het_num, hom_alt_num = get_var_stat(
#         genotype_matrix_chunk)
#     return mis, maf, het, hom_ref_num, het_num, hom_alt_num


def get_var_stat_num_parallel(genotype_matrix, chunk_size=1000, n_jobs=8):
    # Split the genotype_matrix into chunks
    mis, maf, het, hom_ref_num, het_num, hom_alt_num = [], [], [], [], [], []

    total_chunks = genotype_matrix.shape[0] // chunk_size + \
        (1 if genotype_matrix.shape[0] % chunk_size != 0 else 0)
    batch_size = n_jobs * 100
    processed_chunks = 0

    for start_chunk in range(0, total_chunks, batch_size):
        end_chunk = min(start_chunk + batch_size, total_chunks)

        args_dict = {}
        for k in range(start_chunk, end_chunk):
            k = k * chunk_size
            k_end = min(k + chunk_size, genotype_matrix.shape[0])
            genotype_matrix_chunk = genotype_matrix[k:k_end]
            args_dict[k] = (genotype_matrix_chunk,)

        mlt_dict = multiprocess_running(
            get_var_stat, args_dict, n_jobs, silence=True)

        for k in range(start_chunk, end_chunk):
            k = k * chunk_size
            mis_chunk, maf_chunk, het_chunk, hom_ref_num_chunk, het_num_chunk, hom_alt_num_chunk = mlt_dict[
                k]['output']
            mis.extend(mis_chunk)
            maf.extend(maf_chunk)
            het.extend(het_chunk)
            hom_ref_num.extend(hom_ref_num_chunk)
            het_num.extend(het_num_chunk)
            hom_alt_num.extend(hom_alt_num_chunk)

        log_print("processed %d/%d chunks, %.2f%%" % (processed_chunks,
                                                      total_chunks, processed_chunks/total_chunks*100))
        
        processed_chunks += batch_size

    return np.array(mis), np.array(maf), np.array(het), np.array(hom_ref_num), np.array(het_num), np.array(hom_alt_num)


def build_var_stat_table(input_vcf_file, ref_genome_file, output_h5_file):
    chr_len_dict = get_chromosome_info(ref_genome_file)

    for chr_id in chr_len_dict:
        log_print(f'Processing {chr_id}')
        chr_len = chr_len_dict[chr_id]

        genotype_matrix = []
        vcf = VCF(input_vcf_file)
        records = []
        for record in vcf(f'{chr_id}:1-{chr_len}'):
            # chr_var_df = chr_var_df.append({
            #     'ID': record.ID,
            #     'CHROM': record.CHROM,
            #     'POS': record.POS,
            #     'REF': record.REF,
            #     'ALT': record.ALT[0],
            # }, ignore_index=True)
            records.append({
                'ID': record.ID,
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': record.ALT[0],
            })

            variant_genotypes = [
                gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
            genotype_matrix.append(variant_genotypes)

        chr_var_df = pd.DataFrame(records)

        genotype_matrix = np.array(genotype_matrix)

        mis, maf, het, hom_ref_num, het_num, hom_alt_num = get_var_stat_num_parallel(
            genotype_matrix, chunk_size=1000, n_jobs=8)

        chr_var_df['MISSF'] = mis
        chr_var_df['MAF'] = maf
        chr_var_df['HETF'] = het
        chr_var_df['HOM_REF'] = hom_ref_num
        chr_var_df['HET'] = het_num
        chr_var_df['HOM_ALT'] = hom_alt_num

        chr_var_df.to_hdf(output_h5_file, key=chr_id, mode='a')

        vcf.close()
        del genotype_matrix
        del chr_var_df
        gc.collect()

    return output_h5_file


def get_chr_list_from_var_stat_h5(var_stat_h5_file):
    with pd.HDFStore(var_stat_h5_file) as store:
        chr_list = store.keys()
    chr_list = [chr_id.split('/')[1] for chr_id in chr_list]
    return chr_list


if __name__ == '__main__':
    from yxutil import read_list_file

    varID_list_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep"
    varID_list = read_list_file(varID_list_file)

    input_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.vcf.gz"
    output_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/landraces/Sbv5.1.landraces.snp.win10000.maf0.10.miss0.50.rq0.50.ld.nr.rep.vcf"

    output_vcf_file = extract_subvcf_by_varIDs(
        input_vcf_file, varID_list, output_vcf_file)

    genotype_matrix = get_genotype_matrix_from_vcf(output_vcf_file)

    test_vcf_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/test/test.vcf.gz"
    sample_list = get_sample_list_from_vcf(test_vcf_file)
    varID_list = get_varID_list_from_vcf(
        test_vcf_file, chr_id='Chr01', start=1, end=1000)

    test_vcf_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples.vcf.gz'
    ref_genome_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa'
    var_stat_h5_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.georef/population_structure/reseq_africa_landraces/target_samples.var_stat.h5'

    build_var_stat_table(test_vcf_file, ref_genome_file, var_stat_h5_file)
    chr01_var_df = pd.read_hdf(var_stat_h5_file, key='Chr01')
