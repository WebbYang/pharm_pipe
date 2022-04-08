import subprocess
import re
import pandas as pd
import yaml


def call_cmd(cmd):
    '''
    Call shell command in python
    '''
    res = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE).communicate()
    return res

def run_merge_results(out_dir, pid, gender):
    new_df_aldy = pd.read_csv(f'{out_dir}/aldy_{pid}_summary.csv',index_col=1)
    new_df_aldy = new_df_aldy[new_df_aldy.columns[1:]]
    new_df_aldy['CYP2D6'] = new_df_aldy['CYP2D6'].apply(lambda x:x.replace('.ALDY',''))
    df_stargazer = pd.read_csv(f'{out_dir}/stargazer_{pid}_summary.csv',index_col=1)
    df_stargazer = df_stargazer[df_stargazer.columns[1:]]
    df_cyrius = pd.read_table(f'{out_dir}/cyrius_{pid}.tsv', index_col=0)
    df_cyrius.columns = ['CYP2D6','Filter']
    df_cyrius.index = [item.split('.')[0] for item in df_cyrius.index]
    df_astrolabe = pd.read_csv(f'{out_dir}/astrolabe_{pid}_summary.csv', index_col=0)
    
    df_dic = {'star': df_stargazer,'aldy': new_df_aldy,'astrolabe': df_astrolabe,'cyrius': df_cyrius}
    gene_list = list(set(df_stargazer.columns)-set(['CYP4F2', 'IFNL3', 'SLCO1B1']))
    idx = df_stargazer.index
    res_df_all = saveall_pipe(df_dic, gene_list, idx)

    if gender=='M':
        res_df_all['G6PD']['decide'] = res_df_all['G6PD']['decide'].apply(lambda x:x.split('/')[0]+' (male)')

    return res_df_all

def to_yaml(out_dir, pid, res_df_all):

    to_write = {}
    for gene, df in res_df_all.items():       
        genotypes = []
        genotypes.append(df['decide'].loc[pid])
        to_write[gene] = genotypes

    with open(f'{out_dir}/{pid}_starcall.yaml', 'w') as f:
        yaml.dump(to_write, f)   

def reverse_nums(nums):
    '''
    reverse if star allele number is not in order 
    '''
    left, right = nums
    if int(left)>int(right):
        return f'*{right}+*{left}'
    return f'*{left}+*{right}'

def check_rev(word):
    if word=='CannotCall' or word=='None':
        return word
    pattern = re.compile(r'\*(\d+)')
    word_split = word.split('/')
    nums = []
    for item in word_split:
        nums.append(pattern.findall(item))
    left, right = nums    
    left_len, right_len = len(left), len(right)

    if len(left)==1 and len(right)==1:
        if int(left[0])>int(right[0]):
            if '+' not in word_split[0] and '+' not in word_split[1]:
                return f'*{right[0]}/*{left[0]}'
    else:
        if left_len==2:
            left = reverse_nums(left)
        elif left_len==1:
            if '+' not in word_split[0]:
                left = f'*{left[0]}'
            else:
                left = word_split[0]
        else:
            left = word_split[0]
        if right_len==2:
            right = reverse_nums(right)
        elif right_len==1:
            if '+' not in word_split[1]:
                right = f'*{right[0]}'
            else:
                right = word_split[1]
        else:
            right = word_split[1]
        return f'{left}/{right}'
    return word

def equal(pair):
    a, b = pair
    if a=='CannotCall' or b=='CannotCall':
        return False
    return a==b

from itertools import combinations

def concordance(elements, methods, ret_method=False):
    '''
    elements: list of each row vals
    '''
    comb = combinations([(e,m) for e,m in zip(elements, methods)], 2)
    method_record = []
    concord_val = []
    for item in comb:
        es, ms = [i[0] for i in item], [i[1] for i in item]
        if equal(es):
            if es[0] not in concord_val:
                concord_val.append(es[0])
            for m in ms:
                if m not in method_record:
                    method_record.append(m)       
        
    if ret_method:
        em_map = {m:e for e,m in zip(elements, methods)}
        if 'aldy' in method_record:
            # if aldy is in concordance, first priority
            return em_map['aldy']
        try:
            return em_map[method_record[0]]
        except:
            #input(method_record)
            if 'aldy' in methods:
                if em_map['aldy']!='CannotCall':
                    return em_map['aldy']
                else:
                    em_map['star']
            return em_map['star']
    
    if len(concord_val)>1:
        return int(len(method_record)/2)
    return len(method_record)
    

def concordance_pipe(df_dic, gene, index):    
    dic_gene = {f'{gene}_{method}':df_dic[method][gene].tolist() for method in ['star','aldy','astrolabe','cyrius'] 
                 if gene in df_dic[method].columns} 
    df_new = pd.DataFrame(dic_gene, index=index)
    df_new.fillna('CannotCall', inplace=True)

    for col in df_new.columns:
        df_new[col] = df_new[col].apply(check_rev)
        
    methods = [k.split('_')[1] for k in dic_gene.keys()] 
    df_new['concordance'] = df_new.apply(concordance, args=(methods,), axis=1)
    df_new['decide'] = df_new.apply(concordance, args=(methods, True), axis=1)
    return df_new

def saveall_pipe(df_dic, gene_list, idx):
    res_df_all = {}
    for gene in gene_list:
        df = concordance_pipe(df_dic, gene, idx)
        res_df_all[gene] = df[['decide']]
    return res_df_all