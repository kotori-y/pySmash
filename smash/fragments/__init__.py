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
from rdkit import Chem

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
        self.nBits = nBits if folded else None
        self.folded = folded
        self.maxRadius = maxRadius
        self.minRadius = minRadius
        self.maxFragment = maxFragment
        self.nJobs = nJobs if nJobs >= 1 else None

    def GetCircularFragmentLib(self, mols):
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
            Fragments calculated, each element is a list of list
            return from GetCircularFragment() function
        """        
        func = partial(GetCircularFragment,
                       minRadius=self.minRadius,
                       maxRadius=self.maxRadius,
                       nBits=self.nBits,
                       folded=self.folded,
                       maxFragment=self.maxFragment,
                       disposed=True)

        mols = mols if isinstance(mols, Iterable) else (mols, )

        pool = Pool(self.nJobs)
        fragments = pool.map_async(func, mols).get()
        pool.close()
        pool.join()

        return fragments
    
    def _getIdx(self, subArray, Array):
        """Get index of a list to another list
        *internal only*
        """
        uni = set(Array)&set(subArray)
        idx = [Array.index(x) for x in uni]
        return idx
    
    def GetCircularMatrix(self, mols, minNum=None):
        """return circular matrix
    
        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain circular fragments matrix
        minNum : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        pandas.core.frame.DataFrame
            The fragment matrix of molecules
        """

        fragments = self.GetCircularFragmentLib(mols)
        unique = [frag for allFragments in fragments for frag in allFragments[1]]
        
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
         
        return matrix

    def ShowCircularFragment(self, mol, fragmentIndex):
        """Get the SMARTS and SVG image of circular fragment

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The molecule which contain the aim fragment
        fragmentIndex : int
            The index of aim fragment

        Returns
        -------
        smarts : str
            The SMARTS of fragment, which could be used for screening molecules
        svg : str
            The svg string of fragment
        """ 
        fragments = GetCircularFragment(
                mol, minRadius=self.minRadius, maxRadius=self.maxRadius,
                nBits=self.nBits, folded=self.folded, maxFragment=self.maxFragment,
                disposed=False
            )
        center, radius = fragments[fragmentIndex][0]
        smarts, svg = DrawMorganEnv(mol, center, radius)
        return (smarts, svg)
        
        
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


    def GetPathFragmentLib(self, mols):
        """Calculate path-based fragments

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain path-based fragment

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
                       disposed=True)

        pool = Pool(self.nJobs)
        fragments = pool.map_async(func, mols).get()
        pool.close()
        pool.join()

        return fragments
    
    def _getIdx(self, subArray, Array):
        """Get index of a list to another list
        *internal only*
        """
        uni = set(Array) & set(subArray)
        idx = [Array.index(x) for x in uni]
        return idx
    
    def GetPathMatrix(self, mols, minNum=None):
        """return path-based matrix

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain path-based matrix
        minNum : bool, optional
            Whether output with a svg image, by default False

        Returns
        -------
        matrix : pandas.core.frame.DataFrame
            The fragment matrix of molecules
        """
        fragments = self.GetPathFragmentLib(mols)
        
        unique = [bit for info in fragments for bit in info[1]]

        if minNum:
            _all = [bit for info in fragments for bit in info[0]]
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
        for idx,arr in zip(st, matrix):
            arr[idx]=1
        matrix = pd.DataFrame(matrix, columns=unique)
         
        return matrix

    def ShowPathFragment(self, mol, fragmentIndex):
        """Get the SMARTS and SVG image of Path fragment

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The molecule which contain the aim fragment
        fragmentIndex : int
            The index of aim fragment

        Returns
        -------
        smarts : str
            The SMARTS of fragment, which could be used for screening molecules
        svg : str
            The svg string of fragment
        """        
        fragments = GetPathFragment(
                mol, minPath=self.minPath, maxPath=self.maxPath,
                nBits=self.nBits, folded=self.folded, maxFragment=self.maxFragment,
                disposed=False
            )
        
        bondPath = fragments[fragmentIndex][0]
        try:
            smarts, svg = DrawRDKitEnv(mol, bondPath)
        except Chem.KekulizeException:
            kSmi = Chem.MolToSmiles(mol, kekuleSmiles=True)
            mol = Chem.MolFromSmarts(kSmi)
            smarts, svg = DrawRDKitEnv(mol, bondPath)
        return smarts, svg
    

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
        self.nJobs = nJobs if nJobs > 0 else None

    def GetFunctionGroupFragmentLib(self, mols):
        """Calculate function-group fragments

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain function-group fragment

        Returns
        -------
        fragments : list
            Fragments calculated, each element is rerurned from GetFunctionGroupFragment() function
        """
        mols = mols if isinstance(mols, Iterable) else (mols,)
        func = partial(GetFunctionGroupFragment)
        pool = Pool(self.nJobs)
        fgs = pool.map_async(func, mols).get()
        pool.close()
        pool.join()
        return fgs
     
    def _getIdx(self, subArray, Array):
        """
        """
        uni = set(Array) & set(subArray)
        idx = [Array.index(x) for x in uni]
#        print(idx)
        return idx
    
    def GetFunctionGroupsMatrix(self, mols, minNum=None):
        """return function-group matrix

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            The compounds used to obtain function-group matrix

        Returns
        -------
        pandas.core.frame.DataFrame
            the fragment matrix of molecules
        """
        
        fgs = self.GetFunctionGroupFragmentLib(mols)
#        return fgs
#        print(fgs)
    
        unique = [x for fg in fgs for x in fg]
        if minNum:
            c = Counter(unique)
            unique = [x for x in unique if c[x] >= minNum]
        unique = list(set(unique))
        
        num = len(unique)
        matrix = np.zeros((len(mols), num), dtype=np.int8)
        
#        cols =  ['SMARTS'] if not svg else ['SMARTS', 'Substructure']
#        fgs = pd.DataFrame([fg for fg in sum(fgs,[])], columns=cols)
        
#        fgs = pd.DataFrame(fgs)
        # sub = [[x[0] for x in fg] for fg in fgs]
        func = partial(self._getIdx, Array=unique)
        
        pool = Pool(self.nJobs)
        st = pool.map_async(func, fgs).get()
        pool.close()
        pool.join()
#        st = st.apply(lambda x: [dic[i] for i in x])
        for idx,arr in zip(st, matrix):
            arr[idx]=1
        matrix = pd.DataFrame(matrix, columns=unique)

        return matrix

    def ShowFgFragment(self, mol, fragment):
        """Get the SMARTS and SVG image of functional group fragment

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The molecule which contain the aim fragment
        fragment : int
            The fragment in SMARTS format

        Returns
        -------
        smarts : str
            The SMARTS of fragment, which could be used for screening molecules
        svg : str
            The svg string of fragment
        """ 
        svg = DrawFgEnv(mol, fragment)
        smarts = fragment
        return smarts, svg

if '__main__' == __name__:
    freeze_support()
    from rdkit import Chem
    import time
    
    data = pd.read_csv(r'..\..\tests\Carc\Carc.txt', sep='\t')
    # data = data.sample(frac=0.2)
    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    # mols = [Chem.MolFromSmiles(smi) for smi in smis]

    # start = time.clock()
    # circular = Circular(folded=False, maxRadius=7,
    #                     minRadius=1, maxFragment=True,
    #                     nJobs=20)
    # # cirMatrix = circular.GetCircularMatrix(mols, minNum=5)
    # svg = circular.ShowCircularFragment(mols[1115], 2380084179)
    # print(svg)
    # end = time.clock()
    # print(end - start)
    # circular_matrix = circular.GetCircularMatrix(mols, svg=True)

    # path = Path(minPath=1, maxPath=4)
    # # pathMatrix = path.GetPathMatrix(mols)
    # # print(pathMatrix)
    # svg = path.ShowPathFragment(mols[29], 1003675631)
    # print(svg)
    # print(path.substructure)
    # path_matrix = path.GetPathMatrix()
    # print(pathMatrix.shape)

    funcgroup = FunctionGroup()
    print(funcgroup.GetFunctionGroupsMatrix(mols, minNum=5))
    print('Done!!')