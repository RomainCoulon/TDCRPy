# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 08:57:37 2023

@author: jialin.hu
"""

import tdcrpy.TDCR_model_lib as lib

test = lib.readPenNuc2('Co-60')
print(test[-2])