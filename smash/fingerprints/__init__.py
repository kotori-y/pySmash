# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:12:50 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from functools import partial
from collections.abc import Iterable
from multiprocessing import Pool
import numpy as np
import pandas as pd

from .daylight import CalculateDaylight, CalculateSparseDaylight, SmashMolWithDaylight
from .morgan import CalculateMorgan, CalculateSparseMorgan, SmashMolWithMorgan


class Morgan(object):
    
    def __init__(self, mols, 
                 radius=2, nBits=1024, 
                 sparse=False, n_jobs=1):
        """
        Init

        Parameters
        ----------
        mols : a list of rdkit.Chem.rdchem.Mol
            the aim molecule library.
        radius : int, optional
            the radius of circular fingerprints. The default is 2.
            This param would be ignored if sparse set as True
        nBits : int, optional
            the number of bit of morgan. The default is 1024.
            this param would be ignored, if the sparse set as True.
        sparse : bool, optional
            which generate fragment based on sparse fingerprint. 
            The default is False.
        n_jobs : int, optional
            The number of CPUs to use to do the computation
            The default is 1

        Returns
        -------
        None.

        """
        self.mols = mols if isinstance(mols, Iterable) else (mols,)
        self.radius = radius
        self.nBits = nBits if not sparse else 'sparse'
        self.sparse = sparse
        self.n_jobs = n_jobs if n_jobs >=1 else None
        
    def GetMorganFingerAndBitInfo(self):
        """
        Calculate morgan fingerprint and return the info of NonzeroElements

        Returns
        -------
        bitInfo : list of dict
            bit info.

        """
        if self.sparse:
            func = partial(CalculateSparseMorgan, radius=self.radius)
        else:
            func = partial(CalculateMorgan, 
                           radius=self.radius, nBits=self.nBits)
        
        pool = Pool(self.n_jobs)
        bitInfo = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()
        
        self.bitInfo = bitInfo
        return bitInfo
    
    def GetMorganMatrix(self):
        
        func = partial(SmashMolWithMorgan, 
                       radius=self.radius, nBits=self.nBits,
                       sparse=self.sparse)
                       
        
        pool = Pool(self.n_jobs)
        substructures = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()
        
        unique = list(set(sum(substructures, [])))
        dic = dict(zip(unique, range(len(unique))))
        
        num = len(unique)
        
        matrix = np.zeros((len(self.mols), num))
        for idx, arr in enumerate(matrix):
            substructure = substructures[idx]
            for item in substructure:
                arr[dic[item]] = 1
            
        matrix = pd.DataFrame(matrix, columns=unique)
        return matrix
        
        
    
class Daylight(object):
    
    def __init__(self, mols, 
                 minPath=1, maxPath=7, 
                 nBits=1024, sparse=False, 
                 n_jobs=1):
        """
        Init

        Parameters
        ----------
        mols : a list of rdkit.Chem.rdchem.Mol
            the aim molecule library.
        minPath : int, optional
            minimum number of bonds to include in the subgraphs. The default is 1.
        maxPath : int, optional
            maximum number of bonds to include in the subgraphs. The default is 1.
        nBits : int, optional
            the number of bit of daylight. The default is 2048.
            this param would be ignored, if the sparse set as True.
        sparse : bool, optional
            which generate fragment based on sparse fingerprint. 
            The default is True.
        n_jobs : int, optional
            The number of CPUs to use to do the computation
            The default is 1

        Returns
        -------
        None.

        """
        self.mols = mols if isinstance(mols, Iterable) else (mols,)
        self.minPath = minPath
        self.maxPath = maxPath
        self.nBits = nBits if not sparse else 'sparse'
        self.sparse = sparse
        self.n_jobs = n_jobs if n_jobs >=1 else None
        
    def GetDaylightFingerAndBitInfo(self):
        """
        Calculate morgan fingerprint and return the info of NonzeroElements

        Returns
        -------
        fpInfo : list of tuple
            in each tuple, the firsrt element is fingerprint,
            the second is bit info.

        """
        if self.sparse:
            func = partial(CalculateSparseDaylight, 
                           minPath=self.minPath, maxPath=self.maxPath)
        else:
            func = partial(CalculateDaylight,
                           minPath=self.minPath, maxPath=self.maxPath,
                           nBits=self.nBits)
        
        pool = Pool(self.n_jobs)
        bitInfo = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()
        
        self.bitInfo = bitInfo
        return bitInfo
    
    def GetDaylightMatrix(self):
        
        func = partial(SmashMolWithDaylight, 
                       minPath=self.minPath, maxPath=self.maxPath,
                       nBits=self.nBits, sparse=self.sparse)
                       
        
        pool = Pool(self.n_jobs)
        substructures = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()
        
        unique = list(set(sum(substructures, [])))
        dic = dict(zip(unique, range(len(unique))))
        
        num = len(unique)
        
        matrix = np.zeros((len(self.mols), num))
        for idx, arr in enumerate(matrix):
            substructure = substructures[idx]
            for item in substructure:
                arr[dic[item]] = 1
            
        matrix = pd.DataFrame(matrix, columns=unique)
        return matrix

        
        
if '__main__' == __name__:
    from rdkit import Chem
    
    smis = [
            'C1=CC=CC(C(Br)C)=C1',
            'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
            'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
            'C1=NC(CCN)=CN1',
            'C1CCCC(CCO)C1',
            'C1=CC=C2N=C(O)C=CC2=C1',
            'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C',
            'C1=C2N=CC=NC2=C2N=CNC2=C1',
            'C1=C(O)C=CC(O)=C1',
            ]
    
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    
    daylight = Daylight(mols, sparse=True)
    daylight_matrix = daylight.GetDaylightMatrix()
    
    morgan = Morgan(mols, sparse=True)
    morgan_matrix = morgan.GetMorganMatrix()
    
    