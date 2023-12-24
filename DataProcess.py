import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import argparse
from utils.FCGR import CGR


parser = argparse.ArgumentParser()


parser.add_argument('-bacteria', type=str, default='E.coli')
parser.add_argument('-data_dir', type=str, default='/data/HWK/DeepGene/data')
parser.add_argument('-check_data', type=bool, default=False)
parser.add_argument('-drug_name_list', type=list, default=['AMP','AMX','AMC','TZP','CXM','CET','TBM','TMP'])
parser.add_argument('-mode_list', type=list, default=['ToN'])

opt = parser.parse_args()


def check_data(raw_data_dir, download_url_file):
    fa_list = os.listdir(raw_data_dir)
    URL = pd.read_csv(download_url_file, sep=',', encoding='utf-8')
    acce2title = {URL['sample_accession'][i]: URL['sample_title'][i] for i in range(len(URL))}
    errs = []
    for f_name in fa_list:
        with open(raw_data_dir+'/'+f_name, 'r', encoding='utf-8') as f:
            s = re.search(r'SAMEA.*\.', f.readline()).group()
            if acce2title[s[:-1]] != f_name[:-3]:
                errs.append(f_name)

    if len(errs) != 0:
        for e in errs:
            print(e)
    else:
        print("No errors!")


def data_encode(data, label, drug_name, save_dir, mode='raw'):
    """
    mode:
        raw:直接编码
        ToN:将未突变的转化N再进行编码
        combine:将snp序列与参考序列结合构成新序列
    """

    pos = data['POS'].values
    assert (mode in ['raw', 'ToN', 'combine']), "不是合法的mode"

    if mode == 'raw':
        seq_base = {s: i for i, s in enumerate('ATCG')}
    if mode == 'ToN':
        seq_base = {s: i for i, s in enumerate('ATCGN')}
        for name in data.columns[3:]:
            data[name][data[name] == data['REF']] = 'N'
    if mode == 'combine':
        seq_base = {}
        for name in data.columns[3:]:
            data[name] = data['REF']+data[name]
            new_keys = data[name].value_counts().keys()
            for k in new_keys:
                if k not in seq_base:
                    seq_base[k] = len(seq_base)
    # print(seq_base)
    cgr_data = data[data.columns[3:]]
    FCGR_save_dir = os.path.join(save_dir, 'FCGR')
    if not os.path.exists(FCGR_save_dir):
        os.mkdir(FCGR_save_dir)
    cgr_data.to_csv(os.path.join(FCGR_save_dir, f'{opt.bacteria}_{mode}_{drug_name}_FCGR_input.csv'), index=False, sep=',', encoding='utf-8')
    numerical_data = cgr_data.T
    for key, value in seq_base.items():
        numerical_data[numerical_data == key] = value
    numerical_data = numerical_data.values
    Label_Encoding_save_dir = os.path.join(save_dir, 'Label_Encoding')
    if not os.path.exists(Label_Encoding_save_dir):
        os.mkdir(Label_Encoding_save_dir)
    np.savez(os.path.join(Label_Encoding_save_dir, f'{opt.bacteria}_{mode}_Label_Encoding_{drug_name}.npz'), X=numerical_data, Y=label, POS=pos)
    return cgr_data, seq_base


def FCGR_encode(data, label, drug_name, save_dir, mode, seq_base, res=100):
    FCGR_save_dir = os.path.join(save_dir, 'FCGR')
    if data is None:
        data = pd.read_csv(os.path.join(FCGR_save_dir, f'{opt.bacteria}_{mode}_FCGR_input.csv'), sep=',', encoding='utf-8')
    encoded_data = []
    for i, name in enumerate(data.columns):
        print(f"正在处理第{i+1}--{name}个...")
        result = CGR(data=data[name].values, seq_base=seq_base, res=res)
        if i == 0:
            plt.imsave(os.path.join(FCGR_save_dir, f'{opt.bacteria}_{mode}_FCGR{i+1}.png'), result['FCGR'], cmap='gray')
        encoded_data.append(result['FCGR'])

    np.savez(os.path.join(FCGR_save_dir, f'{opt.bacteria}_{mode}_FCGR_{drug_name}.npz'), X=np.array(encoded_data), Y=label)


def Get_label(label_file, sample_list, drug_name):
    pheno2digit = {'R': 1, 'I': 0, 'S': 0}
    phenotypic = pd.read_csv(label_file, sep=',', encoding='utf-8')
    drug = phenotypic[['Isolate', drug_name]]
    drug = drug.dropna(axis=0, how='any').reset_index(drop=True)
    id2label = {drug['Isolate'][i]: pheno2digit[drug[drug_name][i]] for i in range(len(drug))}
    labels = []
    err_sample = []
    for sample in sample_list:
        if sample in id2label:
            labels.append(id2label[sample])
        else:
            err_sample.append(sample)
    return np.array(labels), err_sample


def main():
    data_dir = os.path.join(opt.data_dir, opt.bacteria)
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    raw_data_dir = os.path.join(data_dir, 'raw_data')
    if not os.path.exists(raw_data_dir):
        os.mkdir(raw_data_dir)
    # print(raw_data_dir+'f.txt')
    download_url_file = os.path.join(data_dir, 'download_url.csv')
    # print(download_url_file)

    if opt.check_data:
        check_data(raw_data_dir, download_url_file)

    core_snp_file = os.path.join(data_dir, 'core/core.tab')
    # print(core_snp_file)
    print('Reading raw data ......')
    snp_data = pd.read_csv(core_snp_file, sep='\t', encoding='utf-8')
    print(f"raw_data shape:{snp_data.shape}")

    sample_list = list(snp_data.columns[3:])
    phenotypic_file = os.path.join(data_dir, 'phenotypic.csv')

    preprocess_dir = os.path.join(data_dir, 'preprocessed')
    if not os.path.exists(preprocess_dir):
        os.mkdir(preprocess_dir)

    for drug_name in opt.drug_name_list:
        labels, err_sample = Get_label(phenotypic_file, sample_list, drug_name)
        clear_data = snp_data.drop(columns=err_sample).reset_index(drop=True)

        print(f"raw_data shape:{snp_data.shape} clear_data shape:{clear_data.shape} labels shape:{labels.shape}")
        for mode in opt.mode_list:
            print(f"Encoding for {drug_name} in {mode} mode......")
            # data_encode(clear_data.copy(), labels, drug_name, preprocess_dir, mode=mode)
            cgr_data, seq_base = data_encode(clear_data.copy(), labels, drug_name, preprocess_dir, mode=mode)
            print(f"FCGR Encoding fro {drug_name} in {mode} mode......")
            FCGR_encode(cgr_data, labels, drug_name, preprocess_dir, mode=mode, seq_base=list(seq_base.keys()))
    print("Done!")


if __name__ == '__main__':
    main()
