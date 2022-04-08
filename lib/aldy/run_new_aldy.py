
import subprocess
import argparse
import pathlib
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-b", "--bamdir", required=True)
parser.add_argument("-o", "--outdir", required=True)

args = parser.parse_args()
pid = args.sample

genes = ['CFTR', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'DPYD', 'G6PD', 'IFNL3','NUDT15', 'SLCO1B1', 'TPMT', 'UGT1A1','VKORC1']

sh_path = pathlib.Path(__file__).parent.resolve()

for gene in genes:
    cmd = f'bash {sh_path}/aldy_run.sh {gene} {pid} {args.bamdir} {args.outdir}/aldy'
    result = subprocess.run(cmd, shell=True, universal_newlines=True, check=True)
    print(result.stdout)

def parse_aldy_diplotype(file):
    fh = open(file, 'r')
    for line in fh:
        if line.startswith('#Solution 1'): #Solution 2
            break
    line = fh.readline()
    fh.close()
    try:
        return line.split('\t')[3]
    except:
        return 'NaN'

aldy_summary = {'sample_id':[]}
for gene in genes:
    aldy_summary[gene] = []
    
aldy_summary['sample_id'].append(pid)
for gene in genes:
    file = f'{args.outdir}/aldy/{gene}/{pid}.{gene}.aldy'
    #print(file)
    aldy_summary[gene].append(parse_aldy_diplotype(file))
    
df = pd.DataFrame(aldy_summary)
df.to_csv(f'{args.outdir}/aldy_{pid}_summary.csv')