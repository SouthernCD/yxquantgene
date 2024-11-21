# yxquantgene
Xu Yuxing's personal quantitative genomic tools

## Installation
```bash
pip install yxquantgene
```

## Usage

### 1. Read VCF file and get variant statistics

```python
from yxquantgene import build_var_stat_table

vcf_file = 'path/to/vcf_file' # vcf file should be bgzip compressed and indexed by tabix
ref_genome_file = 'path/to/ref_genome_file' 
stat_h5_file = 'path/to/output_stat_h5_file'

build_var_stat_table(vcf_file, ref_genome_file, stat_h5_file)
```

### 2. Filter variants by statistics

```python
import pandas as pd
from yxquantgene import extract_subvcf_by_varIDs

chr_id = 'Chr01'
max_missing_rate = 0.5
min_maf = 0.01
max_het_rate = 0.1

var_df = pd.read_hdf(stat_h5_file, chr_id)

var_df = var_df[(var_df['MISSF'] <= max_missing_rate)]
var_df = var_df[(var_df['MAF'] >= min_maf)]
var_df = var_df[(var_df['HETF'] <= max_het_rate)]

filtered_varID_list_file = 'path/to/filtered_varID_list_file'
var_df['ID'].to_csv(filtered_varID_list_file, index=False)

filtered_vcf_file = 'path/to/filtered_vcf_file'
extract_subvcf_by_varIDs(input_vcf_file, varID_list, filtered_vcf_file)
```

### 3. Prune variants by LD

You can prune variants by LD and filter low quality variants at the same time.

```python
from yxquantgene import build_LD_db, build_var_stat_table, psa_snp_pruner

input_vcf_file = 'path/to/vcf_file'
ref_genome_file = 'path/to/ref_genome_file'
stat_h5_file = 'path/to/output_stat_h5_file'
snp_ld_dir = 'path/to/snp_ld_dir'
ld_db_win_size = 150000
ld_decay_size = 150000
ld_r2_threshold = 0.5
max_missing_rate = 0.5
min_maf = 0.01
max_het_rate = 0.1
threads = 20
output_prefix = 'path/to/output_prefix'

build_var_stat_table(input_vcf_file, ref_genome_file, stat_h5_file)
build_LD_db(input_vcf_file, stat_h5_file, snp_ld_dir, window_size=ld_db_win_size)
psa_snp_pruner(input_vcf_file, stat_h5_file, ld_db_path, output_prefix, ld_db_win_size=ld_db_win_size, ld_decay_size=ld_decay_size, ld_r2_threshold=ld_r2_threshold, max_missing_rate=max_missing_rate, min_maf=min_maf, max_het_rate=max_het_rate, threads=threads)
```