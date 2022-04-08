import argparse
import os
from lib.utils import *

'''
example: python3 star_call.py --id MMI001 --out out/MMI001 --sex M
'''

parser = argparse.ArgumentParser()
parser.add_argument('--id',help='input subject id', required=True)
parser.add_argument('--out', help='output directory', required=True)
parser.add_argument('--sex', help='M of F', required=True)


args = parser.parse_args()
subject_id = args.id
out_dir = args.out
gender = args.sex

# input file
bam_dir = '/home/dna/webb/low_coverage_test/open_bam'
bam_file = f'{bam_dir}/{subject_id}.realigned.bam'
vcf_dir = '/home/dna/webb/low_coverage_test/open_vcf'
vcf_file = f'{vcf_dir}/{subject_id}.hg19_multianno.vcf'

# output files
bam_file_record = f'{out_dir}/target_bam_file.txt'

# Cyrius
cyrius_dir = "/home/rna/webb/pharm/Cyrius/star_caller.py"
def run_cyrius():
    cmd_list = [f"mkdir -p {out_dir}",
                f"echo {bam_file} > {bam_file_record}",
                f"python3 {cyrius_dir} --manifest {bam_file_record}  --genome 19 \
                                       --prefix cyrius_{subject_id} --outDir {out_dir} --thread 2"]
    for cmd in cmd_list:
        _ = call_cmd(cmd)

    print(f"Call Cyrius to {out_dir}/cyrius_{subject_id} done!")

# Astrolabe
astrolabe_dir = "lib/astrolabe-0.8.7.0"
def run_astrolabe():
    out_file = f"{out_dir}/astrolabe_{subject_id}.txt"
    cmd_list = [f"bash {astrolabe_dir}/sample_run.sh {vcf_file} {bam_file} {out_file}",
                f"python3 {astrolabe_dir}/parse_astrolabe_result.py -i {out_file} -o {out_dir}"]
    for cmd in cmd_list:
        _ = call_cmd(cmd)

# Aldy
aldy_dir = "lib/aldy"
def run_aldy():
    cmd_list = [f"python3 {aldy_dir}/run_new_aldy.py -s {subject_id} -b {bam_file} -o {out_dir}"]
    for cmd in cmd_list:
        _ = call_cmd(cmd)
        
# Stargazer
stargazer_dir = "lib/Stargazer_v1.0.8"
def run_stargazer():
    cmd = f"python3 {stargazer_dir}/run_stargazer.py -s {subject_id} -v {vcf_file} -o {out_dir}"
    _ = call_cmd(cmd)

# merge
def merge_results():
    cmds = [f"mkdir -p {out_dir}/gene_concordance",
           f"mv {out_dir}/*_summary.csv {out_dir}/gene_concordance",
           f"mv {out_dir}/cyrius* {out_dir}/gene_concordance"]
    for cmd in cmds:
        _ = call_cmd(cmd)
    result_df = run_merge_results(f"{out_dir}/gene_concordance", subject_id, gender)
    to_yaml(out_dir, subject_id, result_df)

def main():
    run_cyrius()
    run_astrolabe()
    run_aldy()
    run_stargazer()
    merge_results()

if __name__=="__main__":
    main()

