# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:49:40 2023

@author: jialin.hu
"""

import TDCR_model_lib as tl
#import numpy as np

a,b,c,d,e,f =  tl.readEShape('Fe-59')
print(c,e)
#'''
if len(f[0]) == len(f[0]):
    for i in range(len(f[0])):
        print(tl.incer(f[0][i],e[0][i]))
else:
    print('len pas bon')    
#'''

if len(f) == 2:
    print('### 2 fil')
    if len(f[1]) == len(e[1]):
        for i in range(len(f[1])):
            print(tl.incer(f[1][i],e[1][i]))
    else:
        print('len pas bon')  

