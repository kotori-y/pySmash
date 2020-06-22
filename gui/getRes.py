# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:28:28 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from functools import partial
import multiprocessing as mp
import tkinter as tk
import pandas as pd
import numpy as np
pd.set_option('display.max_colwidth', -1)

from rdkit import Chem
# import openbabel as ob

from smash import Morgan, Daylight, Pvalue, ShowResult

# def obsmitosmile(smi):
#     conv = ob.OBConversion()
#     conv.SetInAndOutFormats("smi", "can")
#     conv.SetOptions("K", conv.OUTOPTIONS)
#     mol = ob.OBMol()
#     conv.ReadString(mol, smi)
#     smile = conv.WriteString(mol)
#     smile = smile.replace('\t\n', '')
#     return smile

def loadmols(smi):
    """
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        with open('./error.log', 'a') as f_obj:
                f_obj.write(smi)
                f_obj.write('/n')
        f_obj.close()
            
    return mol

def CalculatePvalue(n, m, ns, ms, name, subPvalue):
    """
    """
    pValue = sum(Pvalue(n, m, ns, ms))
    subPvalue[name] = pValue


def getFingerprintRes(textPad, data, **kwgrs):

    def add(words):
        textPad['state'] = 'normal'
        textPad.insert(tk.END, words)
        textPad['state'] = 'disable'
    
    ############# Load file #############
    smiles_field = kwgrs.get('smiles_field')
    label_field = kwgrs.get('label_field')
    fingerprint = kwgrs.get('fingerprint')
    minRatio = kwgrs.get('minRatio')
    minNum = kwgrs.get('minNum')
    n_jobs = kwgrs.get('n_jobs')
    aimLabel = kwgrs.get('aimLabel')
    try:
        aimLabel = float(aimLabel)
    except:
        pass
    
    smis = data[smiles_field].values
    labels = data[label_field].values
    add('Successed!\n')
    ############# Load file #############
    

    ############# Trans smiles to mol #############
    add('Trans smiles to mol... ')
    pool = mp.Pool(n_jobs)
    mols = pool.map_async(loadmols, smis).get()
    pool.close()
    pool.join()
    
    mols = pd.Series(mols)
    num = len(mols)
    smis = smis[mols.notna().values]
    labels = labels[mols.notna().values]
    mols = mols[mols.notna().values].to_numpy()
    add('Successed! ({1}/{0})\n'.format(num, len(mols)))
    ############# Trans smiles to mol #############


    ############# Obtain Fingerprint Matrix #############
    add('Obtain Fingerprint Matrix... ')
    if fingerprint == 'ECFP':
        morgan = Morgan(mols, 
                        radius=kwgrs.get('radius'), 
                        minRadius=kwgrs.get('minRadius'),
                        sparse=kwgrs.get('sparse'),
                        nBits=kwgrs.get('nBits'),
                        n_jobs=n_jobs)
        subMatrix = morgan.GetMorganMatrix()
    
    elif fingerprint == 'Daylight':
        daylight = Daylight(mols, 
                          minPath=kwgrs.get('minPath'), 
                          maxPath=kwgrs.get('maxPath'),
                          sparse=kwgrs.get('sparse'),
                          nBits=kwgrs.get('nBits'),
                          n_jobs=n_jobs)
        subMatrix = daylight.GetDaylightMatrix()
        # print(subMatrix)
    add('Successed!\n')
    ############# Obtain Fingerprint Matrix #############
    

    ############# Disposed with run param #############
    add('Disposed with run param... ')
    subMatrix = subMatrix.loc[:, subMatrix.sum(axis=0)>=minNum]
    bo = (subMatrix[labels==aimLabel].sum(axis=0)/subMatrix.sum(axis=0))>=minRatio
    subMatrix = subMatrix.loc[:, bo.values]
    # subMatrix['SMILES'] = smis
    add('Successed!\n')
    ############# Disposed with run param #############


    ############# Calculate p-value #############
    add('Calculate p-Value... ')
    n = len(labels)
    m = sum(labels==aimLabel)
    
    subPvalue = mp.Manager().dict()
    pool = mp.Pool(n_jobs)
    for name, val in subMatrix.iteritems():
        ns = sum(val)
        ms = sum(val==aimLabel)
        pool.apply_async(CalculatePvalue, args=(n, m, ns, ms, name, subPvalue))
    pool.close()
    pool.join()

    subPvalue = pd.DataFrame(dict(subPvalue), index=['Val']).T
    subPvalue = subPvalue[subPvalue['Val']<=0.05]
    subMatrix['SMILES'] = smis
    add('Successed!\n')
    ############# Calculate p-value #############

    return subMatrix, subPvalue, labels


