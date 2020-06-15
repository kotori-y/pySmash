# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 17:03:20 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from multiprocessing import Pool
import numpy as np
from .fingerprints import Morgan, Daylight


# class MorganSmasher(Morgan):
    
#     def __init__(self, mols, **kwgrs):
#         Morgan.__init__(self, mols=mols, **kwgrs)
        
#     def GetMorganMatrix():
#     	pass