
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os


def distr(n, r, plot=False):
    coord = np.zeros((n, 2), dtype='float32')
    for i in range(n):
        coord[i, 0] = r*np.sin(np.pi*((2*i+1)/n))
        coord[i, 1] = r*np.cos(np.pi*((2*i+1)/n))
    if plot:
        plt.scatter(coord[:, 0], coord[:, 1])
        plt.show()
    return coord


def CGR(data=None, seq_base=['A', 'G', 'T', 'C', 'R'], sf=False):
    r = 1.
    base_num = len(seq_base)
    if base_num == 4:
        base_coord = np.array([[1, -1], [-1, -1], [-1, 1], [1, 1]])
        base = {seq_base[i]: base_coord[i] for i in range(base_num)}
    else:
        base_coord = distr(base_num, r)
        base = {seq_base[i]: base_coord[i] for i in range(base_num)}
    # print(base)
    if not sf:
        sf = 1-np.sin(np.pi/base_num)/(np.sin(np.pi/base_num)+np.sin(np.pi*(1/base_num+2*((base_num//4)/base_num))))
    data_len = len(data)
    points = []
    pt = np.zeros(2, dtype='float32')
    for i in range(data_len):
        pt = pt+(base[data[i]]-pt)*sf
        points.append(pt)
    return np.array(points)


def plot_cgr(drug_name, bacteria, mode, topk=50, img_size=224, dpi=50, alpha=10, with_pos=True):
    print(f"loading {drug_name} data...")
    snp = pd.read_csv(f'data/{bacteria}/preprocessed/FCGR/{bacteria}_{mode}_{drug_name}_FCGR_input.csv', sep=',', encoding='utf-8')
    data = np.load(f'/data/HWK/DeepGene/data/{bacteria}/preprocessed/FCGR/{bacteria}_{mode}_FCGR_{drug_name}.npz', allow_pickle=True)
    Y = data['Y'].astype('int32')

    plt.rcParams['savefig.dpi'] = dpi
    figsize = img_size/dpi
    plt.rcParams['figure.figsize'] = (figsize, figsize)
    if with_pos:
        with open(f'results/{bacteria}/Label_Encoding/{drug_name}_POS_weight.json', 'r', encoding='utf-8') as f:
            POS_weight = json.load(f)
        POS_weight = np.array(POS_weight)
        index = np.argsort(-POS_weight)[:topk]
        POS = np.zeros(len(POS_weight))
        POS[index] = 1
    save_dir = f'CGR/{drug_name}_CGR_outputs' if with_pos else f'CGR/raw_{drug_name}_CGR_outputs'
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
        os.mkdir(save_dir+'/0')
        os.mkdir(save_dir+'/1')
    for idx, name in enumerate(snp.columns):
        print(f"{drug_name}-正在处理第{idx+1}个--{name}...")
        points = CGR(snp[name], seq_base=["A", "C", "T", "G", "N"])
        if with_pos:
            no_N = (snp[name] != 'N')
            color = POS*no_N
            s = color*alpha+0.1
        else:
            color = np.zeros(len(points))
            s = 0.1
        plt.scatter(points[:, 0], points[:, 1], c=color, cmap='bwr', s=s)
        plt.axis('off')
        plt.savefig(save_dir+f'/{str(Y[idx])}/{name}.png')
        plt.close()


if __name__ == '__main__':
    # Drug_list = ['AMP', 'AMX', 'AMC', 'TZP', 'CXM', 'CET', 'TBM', 'TMP', 'CIP', 'CTX', 'CTZ', 'GEN']
    Drug_list = ['CTX', 'CTZ']
    Bacteria = 'E.coli'
    Mode = 'ToN'
    print(Drug_list)
    for d in Drug_list:
        # plot_cgr(d, Bacteria, Mode, with_pos=False)
        plot_cgr(d, Bacteria, Mode, with_pos=True)
    print("Done!")