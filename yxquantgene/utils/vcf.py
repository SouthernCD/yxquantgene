from cyvcf2 import VCF, Writer
from yxquantgene.utils.utils import get_chr_len_df
from yxutil import cmd_run, log_print, multiprocess_running, rmdir, have_file
import gc
import numpy as np
import pandas as pd
import time
import random


def read_hdf_with_retry(file_path, key, max_retries=5):
    """
    尝试读取 HDF5 文件，如果失败则重试。

    参数:
    file_path (str): HDF5 文件路径。
    key (str): 要读取的键。
    max_retries (int): 最大重试次数。
    wait_time (int): 每次重试前等待的时间（秒）。

    返回:
    pd.DataFrame: 读取的 DataFrame。
    """
    for attempt in range(max_retries):
        try:
            df = pd.read_hdf(file_path, key=key)
            return df
        except:
            print(f"Attempt {attempt + 1} failed")
            time.sleep(random.uniform(1, 10))
    raise Exception(f"Failed to read HDF5 file after {max_retries} attempts")


def extract_subvcf_by_samples(input_vcf_file, sample_list, output_vcf_file):
    tmp_sample_list_file = 'sample.id.list'
    with open(tmp_sample_list_file, 'w') as f:
        for sample in sample_list:
            f.write(sample + '\n')

    if not output_vcf_file.endswith('.gz'):
        output_vcf_file = output_vcf_file + '.gz'

    cmd_run(
        f"bcftools view --threads 20 -c 1 -O z -o {output_vcf_file} -S sample.id.list {input_vcf_file}")
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
    vcf.close()

    return sample_list


def get_chr_id_list_from_vcf(vcf_file):
    """
    Get chromosome ID list from a VCF file.
    """

    vcf = VCF(vcf_file)
    chr_id_list = vcf.seqnames
    vcf.close()

    return chr_id_list


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


def split_vcf_by_chr(vcf_file, output_dir):
    """
    Split a VCF file by chromosome.
    """

    chr_id_list = get_chr_id_list_from_vcf(vcf_file)

    for chr_id in chr_id_list:
        output_vcf_file = f'{output_dir}/{chr_id}.vcf'
        f, o, e = cmd_run(
            f'bcftools view --threads 20 -O z -o {output_vcf_file}.gz {vcf_file} {chr_id}')
        # print(f,o,e)
        cmd_run(f'tabix -p vcf {output_vcf_file}.gz')


def concat_vcfs(vcf_file_list, output_prefix):
    """
    Concatenate multiple VCF files.
    Will sort the variants by position and remove duplicates only keeping the first SNP (Indel which has the same position will be removed).
    """

    vcf_files = ' '.join(vcf_file_list)
    cmd_run(
        f'bcftools concat -a -d snps -Oz -o {output_prefix}.vcf.gz {vcf_files}')
    cmd_run(
        f'bcftools sort -Oz -o {output_prefix}.sorted.vcf.gz {output_prefix}.vcf.gz')
    cmd_run(f'tabix -p vcf {output_prefix}.sorted.vcf.gz')

    vcf_reader = VCF(f'{output_prefix}.sorted.vcf.gz')
    uniq_vcf_file = f'{output_prefix}.uniq.vcf'

    vcf_writer = Writer(uniq_vcf_file, vcf_reader)

    last_chr = None
    chr_pos_dict = {}
    for record in vcf_reader:
        if last_chr != record.CHROM:
            if len(chr_pos_dict) == 0:
                chr_pos_dict.setdefault(record.POS, []).append(record)
                last_chr = record.CHROM
            else:
                print(f'Processing {last_chr}')
                for pos in sorted(chr_pos_dict.keys()):
                    if len(chr_pos_dict[pos]) == 1:
                        vcf_writer.write_record(chr_pos_dict[pos][0])
                    else:
                        used_rec = None
                        for rec in chr_pos_dict[pos]:
                            if rec.is_snp:
                                used_rec = rec
                                break
                        if used_rec is None:
                            used_rec = chr_pos_dict[pos][0]
                        vcf_writer.write_record(used_rec)

                chr_pos_dict = {}
                chr_pos_dict.setdefault(record.POS, []).append(record)
                last_chr = record.CHROM
        else:
            chr_pos_dict.setdefault(record.POS, []).append(record)

    if len(chr_pos_dict) > 0:
        print(f'Processing {last_chr}')
        for pos in sorted(chr_pos_dict.keys()):
            if len(chr_pos_dict[pos]) == 1:
                vcf_writer.write_record(chr_pos_dict[pos][0])
            else:
                used_rec = None
                for rec in chr_pos_dict[pos]:
                    if rec.is_snp:
                        used_rec = rec
                        break
                if used_rec is None:
                    used_rec = chr_pos_dict[pos][0]
                vcf_writer.write_record(used_rec)

    vcf_writer.close()
    vcf_reader.close()

    cmd_run(f'bgzip {uniq_vcf_file}')
    cmd_run(f'tabix -f -p vcf {uniq_vcf_file}.gz')

    rmdir(f'{output_prefix}.vcf.gz')
    rmdir(f'{output_prefix}.sorted.vcf.gz')
    rmdir(f'{output_prefix}.sorted.vcf.gz.tbi')

    return uniq_vcf_file + '.gz'


def quality_control_vcf(input_vcf_file, output_vcf_file, min_maf=0.05, max_het=0.1, max_miss=0.5):
    if output_vcf_file.endswith('.gz'):
        output_vcf_file = output_vcf_file[:-3]

    vcf_reader = VCF(input_vcf_file)
    vcf_writer = Writer(output_vcf_file, vcf_reader)

    for record in vcf_reader:
        het_num = record.num_het
        hom_ref_num = record.num_hom_ref
        hom_alt_num = record.num_hom_alt
        maf = (het_num + 2*hom_alt_num) / \
            (2*(hom_ref_num + het_num + hom_alt_num))
        maf = min(maf, 1-maf)
        het_freq = het_num / (het_num + hom_ref_num + hom_alt_num)
        het_exp_freq = 2 * maf * (1 - maf)

        if maf < min_maf or het_freq > max_het or record.call_rate < 1 - max_miss:
            continue

        vcf_writer.write_record(record)

    vcf_writer.close()
    vcf_reader.close()

    cmd_run(f'bgzip {output_vcf_file}')
    cmd_run(f'tabix -f -p vcf {output_vcf_file}.gz')

    return output_vcf_file + '.gz'


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


def build_var_stat_table(input_vcf_file, output_h5_file):
    if not have_file(input_vcf_file + ".tbi"):
        cmd_run(f'tabix -p vcf {input_vcf_file}')

    chr_id_list = get_chr_id_list_from_vcf(input_vcf_file)

    for chr_id in chr_id_list:
        log_print(f'Processing {chr_id}')

        genotype_matrix = []
        vcf = VCF(input_vcf_file)
        records = []
        for record in vcf(f"{chr_id}"):
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


def get_chr_list_from_var_stat_h5(var_stat_h5_file, max_retries=5):
    for attempt in range(max_retries):
        try:
            with pd.HDFStore(var_stat_h5_file) as store:
                chr_list = store.keys()
            chr_list = [chr_id.split('/')[1] for chr_id in chr_list]
            return chr_list
        except:
            print(f"Attempt {attempt + 1} failed")
            time.sleep(random.uniform(1, 10))
    raise Exception(f"Failed to read HDF5 file after {max_retries} attempts")


def get_chr_length_dict_from_var_stat_h5(var_stat_h5_file):
    chr_list = get_chr_list_from_var_stat_h5(var_stat_h5_file)
    chr_len_dict = {}
    for chr_id in chr_list:
        df = read_hdf_with_retry(var_stat_h5_file, key=chr_id)
        if len(df) == 0:
            continue
        chr_len = df['POS'].max()
        chr_len_dict[chr_id] = chr_len
    return chr_len_dict


def convert_vcf_to_bimbam(input_vcf_file, output_bimbam_file):
    """
    Convert a VCF file to a BIMBAM file.
    """

    vcf = VCF(input_vcf_file)
    bimbam_file_handle = open(output_bimbam_file, 'w')

    num = 0
    for record in vcf:
        variant_genotypes = [
            gt[0] + gt[1] if gt[0] is not None and gt[1] is not None else -1 for gt in record.genotypes]
        variant_genotypes = [str(x) for x in variant_genotypes]
        bimbam_file_handle.write(
            ",".join([record.ID, record.ALT[0], record.REF] + variant_genotypes) + '\n')
        num += 1

        if num % 10000 == 0:
            print(f'Processed {num} variants.')

    vcf.close()
    bimbam_file_handle.close()

    return output_bimbam_file


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

    build_var_stat_table(test_vcf_file, var_stat_h5_file)
    chr01_var_df = pd.read_hdf(var_stat_h5_file, key='Chr01')
