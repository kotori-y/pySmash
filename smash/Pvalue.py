# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:48:16 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


import scipy as sc
from decimal import Decimal



def Pvalue(n, m, ns, ms):
    """
    get P-value

    Parameters
    ----------
    n : int
        n compounds has m compounds with aspecific activity label.
    m : int
        n compounds has m compounds with aspecific activity label.
    ns : int
        a substructure is found in ns compounds.
    ms : int
        the amount of those compounds with the specific activity label is ms.

    Yields
    ------
    res : float
        P-value.

    """
    for i in range(ms, ns+1):
        numerator = Decimal(sc.math.factorial(float(ns)))
        denominatorA = Decimal(sc.math.factorial(float(i))) * Decimal(sc.math.factorial(float(ns-i)))
        denominatorB = (m/n)**float(i)
        denominatorC = (1 - (m/n))**(ns-i)
        p_value = float(numerator/denominatorA) * denominatorB * denominatorC
        yield p_value
        

     
if '__main__' == __name__:
    p_value = sum(Pvalue(50,10,10,5))





