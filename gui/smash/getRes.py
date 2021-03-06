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


from smash import CircularLearner, PathLearner, FunctionGroupLearner, Pvalue
from itertools import compress
import openbabel as ob
from rdkit import Chem
from functools import partial
import multiprocessing as mp
import tkinter as tk
import pandas as pd
import numpy as np
import time
import bz2
import pickle

try:
    pd.set_option('display.max_colwidth', None)
except:
    pd.set_option('display.max_colwidth', -1)


def obsmitosmile(smi):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = ob.OBMol()
    conv.ReadString(mol, smi)
    smile = conv.WriteString(mol)
    smile = smile.replace('\t\n', '')
    return smile


def loadmols(smi):
    """
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles(obsmitosmile(smi))
        if not mol:
            now = time.strftime('%H:%M:%S %d-%m-%y',
                                time.localtime(time.time()))
            with open('./error.log', 'a') as f_obj:
                f_obj.write('{} {}\n'.format(now, smi))
            f_obj.close()

    return mol


def getFingerprintRes(textPad, data, **kwgrs):

    def add(words):
        textPad['state'] = 'normal'
        textPad.insert(tk.END, words)
        textPad['state'] = 'disable'

    ############# Load file #############
    # print(kwgrs)
    smiles_field = kwgrs.get('smiles_field')
    label_field = kwgrs.get('label_field')
    fingerprint = kwgrs.get('fingerprint')
    minRatio = kwgrs.get('minRatio')
    minNum = kwgrs.get('minNum')
    n_jobs = kwgrs.get('n_jobs')
    aimLabel = kwgrs.get('aimLabel')
    pValue = kwgrs.get('pValue')
    accuracy = kwgrs.get('minAcc')
    Bonferroni = kwgrs.get('Bonferroni')
    try:
        aimLabel = float(aimLabel)
    except:
        pass

    smis = data[smiles_field].values
    labels = data[label_field].values
    add('Successed!\n\n')
    ############# Load file #############

    ############# Trans smiles to mol #############
    add('Trans smiles to mol... ')
    pool = mp.Pool(n_jobs)
    mols = pool.map_async(loadmols, smis).get()
    pool.close()
    pool.join()

    bo = list(map(lambda x: x is not None, mols))
    mols = list(compress(mols, bo))
    smis = list(compress(smis, bo))
    labels = pd.Series((compress(labels, bo)))

    add('Successed! ({0}/{1})\n'.format(len(mols), len(bo)))
    add('{} SMILES Can not be recognized, see error.log\n\n'.format(len(bo)-len(mols)))
    ############# Trans smiles to mol #############

    ############# Obtain Fingerprint Matrix #############
    add('Obtain Significant Fragments... ')
    if fingerprint == 'Circular':
        model = CircularLearner(maxRadius=kwgrs.get('radius'),
                                minRadius=kwgrs.get('minRadius'),
                                folded=kwgrs.get('folded'),
                                maxFragment=True,
                                nJobs=n_jobs)

    elif fingerprint == 'Path':
        model = PathLearner(minPath=kwgrs.get('minPath'),
                             maxPath=kwgrs.get('maxPath'),
                             folded=kwgrs.get('folded'),
                             maxFragment=True,
                             nJobs=n_jobs)

    elif fingerprint == 'Function Group':
        model = FunctionGroupLearner(nJobs=n_jobs)

    subPvalue, subMatrix = model.fit(
        mols, labels,
        aimLabel=aimLabel, minNum=minNum, pCutoff=pValue,
        accCutoff=accuracy, Bonferroni=Bonferroni,
    )

    # print(subMatrix)
    add('Successed!\n\n')
    ############# Obtain Fingerprint Matrix #############

    return model, subMatrix, subPvalue


def predict(modelFile, smis):

    mols = list(map(lambda x: Chem.MolFromSmiles(x), smis))
    model = bz2.BZ2File(modelFile, 'rb')
    model = pickle.load(model)
    y_pred, predMatrix = model.predict(mols)
    return y_pred, predMatrix
