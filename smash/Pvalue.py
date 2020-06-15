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


import pandas as pd
import matplotlib.pyplot as plt



def factorial(num):
    if num == 0 or num == 1:
        return 1
    else:
        return (num * factorial(num - 1))
    
    
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
        part_a = factorial(ns)/(factorial(i) * factorial(ns-i))
        part_b = (m/n)**i
        part_c = (1-m/n)**(ns-i)
        res = part_a * part_b * part_c
        yield res
        

def solve(n, m, ns, ms = 1):
    if 0 < sum(Pvalue(n,m,ns,ms)) <= 0.05:
        return ms
    elif sum(Pvalue(n,m,ns,ms)) == 0:
        return None
    else:
        return solve(n, m, ns, ms+1)
    

def draw():
    X = range(100)
    y_acc = pd.Series([round(x*0.8,1) for x in X])
    y_acc[y_acc < 1] = y_acc[y_acc >= 1].iloc[0]
    
    for ratio in [0.5, 0.7, 0.3]:
        y_p = pd.Series([solve(100/ratio, 100, x) for x in X])
        y_p[y_p.isna()] = y_p[y_p.notna()].iloc[0]
        f,ax = plt.subplots()
        ax.plot(X, y_acc, color='purple',lw=1)
        ax.plot(X, y_p, color='red',lw=1)
        ax.set_xticks(range(0,110,10))
        ax.set_xlim([0,100])
        ax.set_ylim([0,80])
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.tick_params(direction='in', labelsize=8, length=4, colors='black')
        plt.show()
        
if '__main__' == __name__:
    draw()





