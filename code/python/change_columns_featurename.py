# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 18:43:35 2022

@author: wenjin
"""

import pandas as pd
import numpy as np
from collections import Counter
import math
import re
import itertools
import os
from numpy import *
from pandas.core.frame import DataFrame 
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import KFold

number_name=['0.0','0.25','0.5','0.75','1.0']
Vander_distribution_P=['VanDerWaals_distribution_P-'+i for i in number_name]
Vander_distribution_N=['VanDerWaals_distribution_N-'+i for i in number_name]
Vander_distribution_H=['VanDerWaals_distribution_H-'+i for i in number_name] 
Vander_distribution=Vander_distribution_P+Vander_distribution_N+Vander_distribution_H
Vander_composition=['VanDerWaals_composition_'+i for i in ['P','N','H']] 
Vander_transition=['VanDerWaals_transition_'+i for i in ['PN','PH','NH']]
Vander_name=Vander_composition+Vander_transition+Vander_distribution

Hydrophobicity_distribution_P=['Hydrophobicity_distribution_P-'+i for i in number_name]
Hydrophobicity_distribution_N=['Hydrophobicity_distribution_N-'+i for i in number_name]
Hydrophobicity_distribution_H=['Hydrophobicity_distribution_H-'+i for i in number_name] 
Hydrophobicity_distribution=Hydrophobicity_distribution_P+Hydrophobicity_distribution_N+Hydrophobicity_distribution_H
Hydrophobicity_composition=['Hydrophobicity_composition_'+i for i in ['P','N','H']] 
Hydrophobicity_transition=['Hydrophobicity_transition_'+i for i in ['PN','PH','NH']]
Hydrophobicity_name=Hydrophobicity_composition+Hydrophobicity_transition+Hydrophobicity_distribution 

polarity_distribution_P=['Polarity_distribution_P-'+i for i in number_name]
polarity_distribution_N=['Polarity_distribution_N-'+i for i in number_name]
polarity_distribution_H=['Polarity_distribution_H-'+i for i in number_name] 
polarity_distribution=polarity_distribution_P+polarity_distribution_N+polarity_distribution_H
polarity_composition=['Polarity_composition_'+i for i in ['P','N','H']] 
polarity_transition=['Polarity_transition_'+i for i in ['PN','PH','NH']]
polarity_name=polarity_composition+polarity_transition+polarity_distribution 

Polarizability_distribution_P=['Polarizability_distribution_P-'+i for i in number_name]
Polarizability_distribution_N=['Polarizability_distribution_N-'+i for i in number_name]
Polarizability_distribution_H=['Polarizability_distribution_H-'+i for i in number_name] 
Polarizability_distribution=Polarizability_distribution_P+Polarizability_distribution_N+Polarizability_distribution_H
Polarizability_composition=['Polarizability_composition_'+i for i in ['P','N','H']] 
Polarizability_transition=['Polarizability_transition_'+i for i in ['PN','PH','NH']]
Polarizability_name=Polarizability_composition+Polarizability_transition+Polarizability_distribution 
name182=Hydrophobicity_name+Vander_name+polarity_name+Polarizability_name
acid="ACDEFGHIKLMNPQRSTVWY"
AAC_name=['AAC_'+ i for i in list(acid)]
dacid=[]
for i in itertools.product(acid,repeat=2):##产生400个变量特征的名字
    dacid.append("".join(i))
dul=dacid+list(acid)
Dul_name=['Dul-AAC_'+i for i in dacid]
nameall=Dul_name+AAC_name+name182

SS3_com=['SS3_composition_'+i for i in ['C','E','H']]+['SS3_transition_'+i for i in ['CC','CE','CH','EC','EE','EH','HC','HE','HH']]
SS3_distri_C=['SS3_distribution_C-'+i for i in number_name]
SS3_distri_E=['SS8_distribution_E-'+i for i in number_name]
SS3_distri_H=['SS8_distribution_H-'+i for i in number_name]
SS3_name=SS3_com+SS3_distri_C+SS3_distri_E+SS3_distri_H

SS8='CSTHGIEB'
doub=[]
for i in itertools.product(SS8,repeat=2):
    doub.append("".join(i))
doub_name=list(SS8)+doub



SS8_com=['SS8_composition_'+i for i in list(SS8)]+['SS8_Dual_'+i for i in doub]
SS8_distri=[]
for i in list(SS8):
    SS8_distri=SS8_distri+['SS8_distribution_'+i+'-'+j for j in number_name]
SS8_name=SS8_com+SS8_distri
other=['ASA','HSEa-u','HSEa-d','CN13','BA_theta','BA_tau','BA_phi','BA_psi','P(3-C)','P(3-E)','P(3-H)']+['P(8-'+ i+')'for i in list(SS8)]
other_name=[]
for i in other:
    other_name=other_name+[i+'-'+j for j in ['1st','2nd','3th']]
structural_name=SS3_name+SS8_name+other_name
featurename=nameall+structural_name
all_feature=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/all_x_matrix.csv').iloc[:,1:]
all_feature.columns=featurename
outpath="F:/zhulin/1datasz1/RNA/outdata"
all_feature.to_csv(os.path.join(outpath,"all_x_matrix_name.csv"))