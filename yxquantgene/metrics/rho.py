import os
from yxutil import mkdir, cmd_run
from yxquantgene.utils.vcf import get_chr_id_list_from_vcf, split_vcf_by_chr

current_script_path = os.path.abspath(__file__)
current_script_dir = os.path.dirname(current_script_path)
fasteprr_script_path = os.path.join(
    current_script_dir, 'Rscripts', 'FASTEPRR.R')


def calulate_rho(input_vcf_file, output_path, window_size=50000):
    """
    Calculate the recombinational parameter rho from a VCF file using the FASTEPRR
    """

    mkdir(output_path)

    clean_vcf_file = os.path.join(output_path, 'input.clean.vcf.gz')
    cmd_str = f'bcftools annotate -x INFO {input_vcf_file} -Oz -o {clean_vcf_file}'
    cmd_run(cmd_str, cwd=output_path)
    cmd_run(f'tabix -p vcf {clean_vcf_file}', cwd=output_path)

    chr_id_list = get_chr_id_list_from_vcf(clean_vcf_file)
    chr_id_list_file = os.path.join(output_path, 'chr_id_list.txt')
    with open(chr_id_list_file, 'w') as f:
        for chr_id in chr_id_list:
            f.write(chr_id + '\n')
    split_vcf_by_chr(clean_vcf_file, output_path)

    cmd_str = f'Rscript {fasteprr_script_path} {output_path} {window_size} {chr_id_list_file}'
    cmd_run(cmd_str, cwd=output_path)


if __name__ == "__main__":
    input_vcf_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/landraces_snp.vcf.gz'
    output_path = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation.2.0/1.envGWAS/population_structure/fasteprr'
