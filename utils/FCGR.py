import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def distr(n, r, plot=False):
    coord = np.zeros((n, 2), dtype='float32')
    for i in range(n):
        coord[i, 0] = r*np.sin(np.pi*((2*i+1)/n))
        coord[i, 1] = r*np.cos(np.pi*((2*i+1)/n))
    if plot:
        plt.scatter(coord[:, 0], coord[:, 1])
        plt.show()
    return coord


def CGR(data=None, seq_base=['A', 'G', 'T', 'C', 'R'], sf=False, res=100):
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
    A = np.zeros((res, res), dtype='float32')
    pt = np.zeros(2, dtype='float32')
    for i in range(data_len):
        pt = pt+(base[data[i]]-pt)*sf
        points.append(pt)
        x = ((pt[0]+r)*res/(2*r)).astype('int')
        y = ((pt[1]+r)*res/(2*r)).astype('int')
        # print(x,y)
        A[x][y] = A[x][y]+1
    return {'FCGR': A, 'cgr_points': points}


Gene=[
    (898515,899645),
    (2655075,2656358),
    (3098558,3099565),
    (4100810,4101430),
    (4172057,4173085),
    (4439875,4441215),
    (4453583,4454578),
    (4466299,4467246),
    (4477307,4478311),
    (4604875,4605663),
    (4627315,4628547)
]

def In_Gene(pos):
    is_in=False
    for x in Gene:
        is_in=is_in|((pos>=x[0])&(pos<=x[1]))
    return is_in

# def CGR2(data=None, seq_base=['A', 'G', 'T', 'C', 'R'], sf=False, res=100):
#     r = 1.
#     base_num = len(seq_base)
#     if base_num == 4:
#         base_coord = np.array([[1, -1], [-1, -1], [-1, 1], [1, 1]])
#         base = {seq_base[i]: base_coord[i] for i in range(base_num)}
#     else:
#         base_coord = distr(base_num, r)
#         base = {seq_base[i]: base_coord[i] for i in range(base_num)}
#     # print(base)
#     if not sf:
#         sf = 1-np.sin(np.pi/base_num)/(np.sin(np.pi/base_num)+np.sin(np.pi*(1/base_num+2*((base_num//4)/base_num))))
#     data_len = len(data)
#     points = []
#     A = np.zeros((res, res), dtype='float32')
#     pt = np.zeros(2, dtype='float32')
#     for i in range(data_len):
#         pt = pt+(base[data[i]]-pt)*sf
#         points.append(pt)
#         x = ((pt[0]+r)*res/(2*r)).astype('int')
#         y = ((pt[1]+r)*res/(2*r)).astype('int')
#         # print(x,y)
#         if In_Gene(POS[i]):
#             # print(x,y)
#             A[x][y] = A[x][y]+10
#         else:
#             A[x][y] = A[x][y]+1
#     return {'FCGR': A, 'cgr_points': points}

