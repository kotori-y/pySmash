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


from collections import Counter
from functools import partial
from collections.abc import Iterable
from multiprocessing import Pool
from multiprocessing import freeze_support
import numpy as np
import pandas as pd

try:
    from .path import *
    from .circular import *
    from .functionGroup import *
except:
    from path import *
    from circular import *
    from functionGroup import *


class Circular:
    """Calculte circular fragments and obtain matrix  
    """    
    def __init__(self, minRadius=1, maxRadius=2, 
                 nBits=1024, folded=False,
                 maxFragment=True, nJobs=1):
        """Initialization

        Parameters
        ----------
        minRadius : int, optional
            The probable minimum radius of circular fragment, by default 1
        maxRadius : int, optional
            The probable maximum radius of circular fragment, by default 2
        nBits : int, optional
            the number of bit of morgan, by default 1014
            this param would be ignored, if the folded set as False.
        folded : bool, optional
            which generate fragment based on unfolded fingerprint, by default True.
        maxFragment : bool, optional
            Whether only return the maximum fragment at a center atom, by default True
        nJobs : int, optional
            The number of CPUs to use to do the computation, by default 1
        """
        # self.mols = mols if isinstance(mols, Iterable) else (mols, )
        self.nBits = nBits if folded else None
        self.folded = folded
        self.maxRadius = maxRadius
        self.minRadius = minRadius
        self.maxFragment = maxFragment
        self.nJobs = nJobs if nJobs >= 1 else None
#        self.substructure = pd.DataFrame()
#        self.matrix = pd.DataFrame()
#        self.fragments = None
    
    def GetCircularFragmentLib(self, mols, svg=False):
        """Calculate circular fragments

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            The compounds used to obtain circular fragments.
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        fragments : list
            Fragments calculated, each element is a tuple 
            return from GetCircularFragment() function
        """        
        func = partial(GetCircularFragment,
                       minRadius=self.minRadius,
                       maxRadius=self.maxRadius,
                       nBits=self.nBits,
                       folded=self.folded,
                       maxFragment=self.maxFragment,
                       svg=svg)

        mols = mols if isinstance(mols, Iterable) else (mols, )

        pool = Pool(self.nJobs)
        fragments = pool.map_async(func, mols).get()
        pool.close()
        pool.join()

        return fragments
    
    def _getIdx(self, subArray, Array):
        """
        """
        uni = set(Array)&set(subArray)
        idx = [Array.index(x) for x in uni]
        return idx
    
    def GetCircularMatrix(self, mols, minNum=None, svg=False):
        """return circular matrix
    
        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain circular fragments matrix
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        pandas.core.frame.DataFrame
            The fragment matrix of molecules
        """

        fragments = self.GetCircularFragmentLib(mols, svg=svg)
        unique = [bit for info in fragments for bit in info[1]]
        
        if minNum:
            _all = [bit for info in fragments for bit in info[0]]
            c = Counter(_all)
            unique = [x for x in unique if c[x] >= minNum]
            
            
        unique = list(set(unique))

        num = len(unique)
        matrix = np.zeros((len(mols), num), dtype=np.int8)
        
        fragments = pd.DataFrame(fragments)
        subArray = fragments[0]
        func = partial(self._getIdx, Array=unique)
        pool = Pool(self.nJobs)
        st = pool.map_async(func, subArray).get()
        pool.close()
        pool.join()
        
        for idx,arr in zip(st, matrix):
            arr[idx]=1
        matrix = pd.DataFrame(matrix, columns=unique)
         
        substructure = {k:v for item in fragments[1].values for k,v in item.items()}
        idx = ['SMARTS'] if not svg else ['SMARTS', 'Substructure']
        substructure = pd.DataFrame(substructure, index=idx).T
        return matrix, substructure


class Path:
    """Calculte path-based fragments and obtain matrix  
    """    
    def __init__(self,
                 minPath=1, maxPath=7,
                 nBits=1024, folded=False,
                 maxFragment=True, nJobs=1):
        
        """Initialization

        Parameters
        ----------
        minPath : int, optional
            The probable minimum length of path-based fragment, by default 1
        maxPath : int, optional
            The probable maximum length of path-based fragment, by default 7
        nBits : int, optional
            the number of bit of morgan, by default 1014
            this param would be ignored, if the folded set as False.
        folded : bool, optional
            which generate fragment based on unfolded fingerprint, by default True.
        maxFragment : bool, optional
            Whether only return the maximum fragment at a center atom, by default True
        nJobs : int, optional
            The number of CPUs to use to do the computation, by default 1
        """
        self.minPath = minPath
        self.maxPath = maxPath
        self.nBits = nBits if folded else None
        self.folded = folded
        self.nJobs = nJobs if nJobs >= 1 else None
        self.maxFragment = maxFragment


    def GetPathFragmentLib(self, mols, svg=False):
        """Calculate path-based fragments

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain path-based fragment
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        fragments : list
            Fragments calculated, each element is returned from GetPathFragment() function
        """
        mols = mols if isinstance(mols, Iterable) else (mols,)
        func = partial(GetPathFragment,
                       minPath=self.minPath,
                       maxPath=self.maxPath,
                       nBits=self.nBits,
                       folded=self.folded,
                       maxFragment=self.maxFragment,
                       svg=svg)

        pool = Pool(self.nJobs)
        bitInfo = pool.map_async(func, mols).get()
        pool.close()
        pool.join()

        return bitInfo
    
    def _getIdx(self, subArray, Array):
        """
        """
        uni = set(Array) & set(subArray)
        idx = [Array.index(x) for x in uni]
#        print(idx)
        return idx
    
    def GetPathMatrix(self, mols, minNum=None, svg=False):
        """return path-based matrix

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain path-based matrix
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        pandas.core.frame.DataFrame
            the fragment matrix of molecules
        """
        fragments = self.GetPathFragmentLib(mols, svg=svg)
        
        unique = [bit for info in fragments for bit in info[1]]

        if minNum:
            _all = [bit for info in fragments for bit in info[0]]
#            _uni = set(unique)&set(_all)
            c = Counter(_all)
            unique = [x for x in unique if c[x] >= minNum]
            
            
        unique = list(set(unique))

        num = len(unique)
        matrix = np.zeros((len(mols), num), dtype=np.int8)
    
        fragments = pd.DataFrame(fragments)
        sub = fragments[0]
        func = partial(self._getIdx, Array=unique)
        
        pool = Pool(self.nJobs)
        st = pool.map_async(func, sub).get()
        pool.close()
        pool.join()
#        st = st.apply(lambda x: [dic[i] for i in x])
        for idx,arr in zip(st, matrix):
            arr[idx]=1
        matrix = pd.DataFrame(matrix, columns=unique)
         
        substructure = {k:v for item in fragments[1].values for k,v in item.items() if k in matrix.columns}
        idx = ['SMARTS'] if not svg else ['SMARTS', 'Substructure']
        substructure = pd.DataFrame(substructure, index=idx).T
        return matrix, substructure


class FunctionGroup:
    """Calculte function-group fragments and obtain matrix  
    """    
    def __init__(self, nJobs=1):
        """Initialization

        Parameters
        ----------
        nJobs : int, optional
            The number of CPUs to use to do the computation, by default 1
        """
        self.nJobs = nJobs
        self.fgs = None
        self.substructure = pd.DataFrame()
        self.matrix = pd.DataFrame()

    def GetFunctionGroupFragmentLib(self, mols, svg=False):
        """Calculate function-group fragments

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain function-group fragment
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        fragments : list
            Fragments calculated, each element is rerurned from GetFunctionGroupFragment() function
        """
        mols = mols if isinstance(mols, Iterable) else (mols,)
        func = partial(GetFunctionGroupFragment, svg=svg)
        pool = Pool(self.nJobs)
        self.fgs = pool.map_async(func, mols).get()
        pool.close()
        pool.join()
        return self.fgs

    def GetFunctionGroupsMatrix(self, mols, svg=False):
        """return function-group matrix

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain function-group matrix
        svg : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        pandas.core.frame.DataFrame
            the fragment matrix of molecules
        """
        if not self.fgs:
            fgs = self.GetFunctionGroupFragmentLib(mols, svg)
        else:
            fgs = self.fgs

        unique = [x[0] for fg in fgs for x in fg]
        unique = list(set(unique))
        dic = dict(zip(unique, range(len(unique))))

        num = len(unique)
        matrix = np.zeros((len(mols), num), dtype=np.int8)

        for idx, arr in enumerate(matrix):
            fg = fgs[idx]
            for item in fg:
                arr[dic[item[0]]] = 1

        self.matrix = pd.DataFrame(matrix, columns=unique)
        cols =  ['SMARTS'] if not svg else ['SMARTS', 'Substructure']
        self.substructure = pd.DataFrame([fg for fg in sum(fgs,[])], columns=cols)
        self.substructure = self.substructure.drop_duplicates(['SMARTS']).set_index('SMARTS', drop=False)
        self.matrix = self.matrix.reindex(self.substructure.SMARTS.values, axis=1)

        return self.matrix


if '__main__' == __name__:
    freeze_support()
    from rdkit import Chem
    import time
    
    data = pd.read_csv(r'..\..\tests\Carc\Carc.txt', sep='\t')
    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    # mols = [Chem.MolFromSmiles(smi) for smi in smis]

    start = time.clock()
    circular = Circular(folded=False, maxRadius=7,
                        minRadius=1, maxFragment=True,
                        nJobs=20)
    fragments = circular.GetCircularMatrix(mols, svg=True)
    end = time.clock()
    print(end - start)
    # circular_matrix = circular.GetCircularMatrix(mols, svg=True)

    # path = Path(minPath=1, maxPath=3)
    # path_frag = path.GetPathMatrix(mols, svg=True)
    # print(path.substructure)
    # path_matrix = path.GetPathMatrix()
    # print(path_matrix.shape)

    # funcgroup = FunctionGroup()
    # print(funcgroup.GetFunctionGroupsMatrix(mols, svg=False))
    # print(funcgroup.substructure)
    print('Done!!')