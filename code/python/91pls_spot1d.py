# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:14:57 2021

@author: wenjin
"""

import pandas as pd

with open("F:/zhulin/1datasz1/RNA/predictRBPspot/outputs/predict_99.spot1d",'r') as f:
    spot_col_names=f.readline().split('\t')
    spot_col_names.pop()
    spot = pd.read_csv(f,delim_whitespace=True,names=spot_col_names)
    spot.drop(axis=0,index=0)
    spot.drop(axis=1,index=0)
