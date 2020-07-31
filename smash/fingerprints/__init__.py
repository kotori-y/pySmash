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


from collections import ChainMap
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


class Circular(object):

    def __init__(self, mols,
                 maxRadius=2, minRadius=1,
                 nBits=1024, folded=True,
                 maxFragment=False, nJobs=1):
        """
        Init

        Parameters
        ----------
        mols : a list of rdkit.Chem.rdchem.Mol
            the aim molecule library.
        maxRadius : int, optional
            the radius of circular fingerprints. The default is 2
        nBits : int, optional
            the number of bit of morgan. The default is 1024
            this param would be ignored, if the folded set as False
        folded : bool, optional
            which generate fragment based on unfoled fingerprint
            The default is True.
        minRadius : int, optional
            minimum radius of environment. The default is 1
        maxFragment : bool, optional
            whether only return the max fragment of each atom, by default False
        nJobs : int, optional
            The number of CPUs to use to do the computation
            The default is 1

        Returns
        -------
        None.

        """
        self.mols = mols if isinstance(mols, Iterable) else (mols, )
        self.nBits = nBits if folded else None
        self.folded = folded
        self.maxRadius = maxRadius
        self.minRadius = minRadius
        self.maxFragment = maxFragment
        self.nJobs = nJobs if nJobs >= 1 else None

    def GetCircularFragmentLib(self):
        """
        Calculate morgan fingerprint and return the info of NonzeroElements

        Returns
        -------
        bitInfo : list of dict
            fragment bit info.

        """
        func = partial(GetCircularFragment,
                       minRadius=self.minRadius,
                       maxRadius=self.maxRadius,
                       nBits=self.nBits,
                       folded=self.folded,
                       maxFragment=self.maxFragment)

        pool = Pool(self.nJobs)
        bitInfo = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()

        self.bitInfo = bitInfo
        return bitInfo

    def GetCircularMatrix(self):
        """return circular matrix

        Returns
        -------
        pandas.core.frame.DataFrame
            the fragment matrix of molecules

        """
        bitInfo = self.GetCircularFragmentLib()
        # pool.close()
        # pool.join()

        unique = [bit for info in bitInfo for bit in info[1]]
        unique = list(set(unique))
        dic = dict(zip(unique, range(len(unique))))

        num = len(unique)
        matrix = np.zeros((len(self.mols), num), dtype=np.int8)

        for idx, arr in enumerate(matrix):
            info = bitInfo[idx]
            for frag in unique:
                if frag in info[0]:
                    arr[dic[frag]] = 1

        # colDict = ChainMap(*[info[1] for info in bitInfo])
        matrix = pd.DataFrame(matrix, columns=unique)
        return matrix


class Path(object):

    def __init__(self, mols,
                 minPath=1, maxPath=7,
                 nBits=1024, folded=False,
                 nJobs=1, maxFragment=True):
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
            this param would be ignored, if the folded set as False.
        folded : bool, optional
            which generate fragment based on sparse fingerprint. 
            The default is True.
        nJobs : int, optional
            The number of CPUs to use to do the computation
            The default is 1

        Returns
        -------
        None.

        """
        self.mols = mols if isinstance(mols, Iterable) else (mols,)
        self.minPath = minPath
        self.maxPath = maxPath
        self.nBits = nBits if folded else None
        self.folded = folded
        self.nJobs = nJobs if nJobs >= 1 else None
        self.maxFragment = maxFragment

    def GetPathFragmentLib(self):
        """return the info of path bit info

        Returns
        -------
        fpInfo : list of tuple
            in each tuple, the firsrt element is fingerprint,
            the second is bit info.

        """
        func = partial(GetPathFragment,
                       minPath=self.minPath,
                       maxPath=self.maxPath,
                       nBits=self.nBits,
                       folded=self.folded,
                       maxFragment=self.maxFragment)

        pool = Pool(self.nJobs)
        bitInfo = pool.map_async(func, self.mols).get()
        pool.close()
        pool.join()

        self.bitInfo = bitInfo
        return bitInfo

    def GetPathMatrix(self):

        bitInfo = self.GetPathFragmentLib()
        # pool.close()
        # pool.join()

        unique = [bit for info in bitInfo for bit in info[1]]
        unique = list(set(unique))
        dic = dict(zip(unique, range(len(unique))))

        num = len(unique)
        matrix = np.zeros((len(self.mols), num), dtype=np.int8)

        for idx, arr in enumerate(matrix):
            info = bitInfo[idx]
            for frag in unique:
                if frag in info[0]:
                    arr[dic[frag]] = 1

        # colDict = ChainMap(*[info[1] for info in bitInfo])
        matrix = pd.DataFrame(matrix, columns=unique)
        return matrix


class FunctionGroup(object):

    def __init__(self, mols, nJobs=1):
        """Init

        Parameters
        ----------
        mols : a list of rdkit.Chem.rdchem.Mol
            the aim molecule library.
        """
        self.mols = mols
        self.nJobs = nJobs
        self.fgs = None

    def GetFunctionGroups(self):

        pool = Pool(self.nJobs)
        self.fgs = pool.map_async(GetFunctionGroupFragment, self.mols).get()
        pool.close()
        pool.join()

        return self.fgs

    def GetFunctionGroupsMatrix(self):

        if not self.fgs:
            fgs = self.GetFunctionGroups()
        else:
            pass

        unique = [x for fg in fgs for x in fg]
        unique = list(set(unique))
        dic = dict(zip(unique, range(len(unique))))

        num = len(unique)
        matrix = np.zeros((len(self.mols), num), dtype=np.int8)

        for idx, arr in enumerate(matrix):
            fg = fgs[idx]
            for item in fg:
                arr[dic[item]] = 1

        matrix = pd.DataFrame(matrix, columns=unique)
        return matrix


if '__main__' == __name__:
    freeze_support()
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
        'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
        'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
        'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C',
        ]

    mols = [Chem.MolFromSmiles(smi) for smi in smis]

    # circular = Circular(mols, folded=False, maxRadius=3,
    #                     minRadius=1, maxFragment=False)
    # circular_matrix = circular.GetCircularMatrix()
    # print(circular_matrix)
    # circular_matrix.insert(0, 'SMILES', smis)
    # circular_matrix.to_csv(r'C:\Users\0720\Desktop\py_work\pySmash\tests\Ames\data0709.csv', index=False)

    path = Path(mols, folded=False, maxFragment=True)
    # path_frag = path.GetPathFragmentLib()
    # print(path_frag[0])
    path_matrix = path.GetPathMatrix()
    print(path_matrix.shape)

    # fg = FunctionGroup(mols)
    # print(fg.GetFunctionGroupsMatrix())
