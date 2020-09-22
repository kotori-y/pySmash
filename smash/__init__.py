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


from decimal import Decimal
import multiprocessing as mp
from collections.abc import Iterable

import os
import pandas as pd
import numpy as np
import scipy as sc
import bz2
import _pickle as cPickle

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor

try:
    from .fragments import Circular, Path, FunctionGroup
except Exception:
    from fragments import Circular, Path, FunctionGroup

try:
    pd.set_option('display.max_colwidth', None)
except:
    pd.set_option('display.max_colwidth', -1)

from scipy.stats import binom


class NotFittedError(Exception):
    pass


def Pvalue(n, m, ns, ms):
    """Get P-value.

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
        use sum() function to get the pvalue.

    """
    p_value = 0
    for i in range(ms, ns+1):
        numerator = Decimal(sc.math.factorial(float(ns)))
        denominatorA = Decimal(sc.math.factorial(
            float(i))) * Decimal(sc.math.factorial(float(ns-i)))
        denominatorB = (m/n)**float(i)
        denominatorC = (1 - (m/n))**(ns-i)
        val = float(numerator/denominatorA) * denominatorB * denominatorC
        p_value += val

    return p_value


class BaseLearner:
    """Base class for all learner in pysamsh.
    """

    def __init__(self):
        """Initialization.
        """
        self.sigFragments = None

    def GetMatrix(self, mols, **kwgrs):
        """Get the matrix of fragment of mol library.

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            compounds, which have aspecific endpoint label, used to obtain significant fragments.
        """
    
    def ShowFragment(self, mol, **kwgrs):
        """Visulaizing fragments.

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            compounds, which have aspecific endpoint label, used to obtain significant fragments.
        """
        pass

    def fit(self, mols,
            labels, aimLabel=1,
            minNum=5, pCutoff=0.05,
            accCutoff=None, Bonferroni=False):
        """Learning from given database which have aspecific endpoint label.

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            Compounds, which have aspecific endpoint label, used to obtain significant fragments
        labels : array-like of shape (len(mols),)
            The target values (class labels in classification)
        aimLabel : any, optional
            The label to be regarded as activity label (class labels in classification), by default 1
        minNum : int, optional
            The minimum frequency a fragment required, by default 5
        pCutoff : float, optional
            The pvalue cutoff, a fragment would be regarded as significant if its pvalue below pCutoff, by default 0.05
        accCutoff : float, optional
            The minimum accraucy lead by a fragment judge, by default None
        Bonferroni : bool, optional
            Whether use Bonferroni method to revised, by default False
        svg : bool, optional
            Whether output with a svg image, by default True

        Returns
        -------
        smash.BaseLearner
            A fitted learner, used predict() method can predict molecules without known label 
        """
        mols = tuple(mols) if isinstance(mols, Iterable) else (mols, )

        matrix = self.GetMatrix(mols, minNum=minNum)

        n = len(labels)
        m = (labels == aimLabel).sum()
        p = m/n

        sigPvalue = pd.DataFrame(
            {"ns": matrix.sum(),
             "ms": matrix[labels == 1].sum()}
        )
        sigPvalue = sigPvalue.apply(
            lambda x: sum(
                binom.pmf(range(x['ms'], x['ns']+1), x['ns'], p)
            ), axis=1
        )

        sigPvalue = sigPvalue[sigPvalue <= pCutoff]

        sigPvalue = pd.DataFrame(sigPvalue, columns=['Pvalue'])
        sigMatrix = matrix.reindex(sigPvalue.index, axis=1)
        sigPvalue['Total'] = sigMatrix.sum(axis=0)
        sigPvalue['Hitted'] = sigMatrix[labels == 1].sum(axis=0)
        sigPvalue['Accuracy'] = sigPvalue.Hitted/sigPvalue.Total
        sigPvalue['Coverage'] = sigPvalue['Hitted']/m

        if accCutoff:
            sigPvalue = sigPvalue[sigPvalue.Accuracy >= accCutoff]
        else:
            pass

        if Bonferroni:
            sigPvalue['Pvalue'] = sigPvalue.Pvalue.values * len(sigPvalue)
            sigPvalue = sigPvalue[sigPvalue.Pvalue <= pCutoff]
        else:
            pass

        substructure = {}
        for idx,vals in sigMatrix.iteritems():
            mol = mols[vals[vals==1].index[0]]
            smarts, svg = self.ShowFragment(mol, idx)
            substructure[idx] = {"SMARTS":smarts, "Substructure":svg}
        substructure = pd.DataFrame(substructure).T
        
        sigPvalue = sigPvalue.merge(substructure, left_index=True, right_index=True)
        

        sigPvalue = sigPvalue.sort_values('Pvalue')
        print(type(sigPvalue))
        sigMatrix = sigMatrix.reindex(sigPvalue.index, axis=1)

        
        return sigPvalue, sigMatrix

    def predict(self, mols):
        """Predict molecules without known label 

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object
            Compounds to be predicted

        Returns
        -------
        y_pred : ndarray of shape (len(mols),)
            The predicted label
        predMatrix : pandas.core.frame.DataFrame
            The predictd matrix 

        Raises
        ------
        NotFittedError
            Judge a learner whether was fitted
        """
        if self.sigFragments is None:
            raise NotFittedError(
                "This instance is not fitted yet. Call 'fit' with appropriate arguments before using this method.")

        predMatrix, _ = self.GetMatrix(mols, svg=True)
        # cols = set(predMatrix.columns) & set(self.sigPvalue.index)
        predMatrix = predMatrix.reindex(
            self.sigFragments.keys(), axis=1).reset_index(drop=True).fillna(0)
        y_pred = ((predMatrix.sum(axis=1) > 0) + 0).values

        return y_pred, predMatrix

    def savePvalue(self, sigPvalue, file='./pvalue.html'):
        """Save pvalue result

        Parameters
        ----------
        file : str, optional
            The path to save result, by default './pvalue.html'
        """
        pd.set_option('colheader_justify', 'center')   # FOR TABLE <th>

        html_string = """<html>
        <head><title>HTML Pandas Dataframe with CSS</title></head>
        <style>
        /* includes alternating gray and white with on-hover color */

        .mystyle {
            font-size: 11pt; 
            font-family: Arial;
            border-collapse: collapse; 
            border: 1px solid silver;

        }

        .mystyle td, th {
            padding: 5px;
        }

        .mystyle tr:nth-child(even) {
            background: #E0E0E0;
        }

        .mystyle tr:hover {
            background: silver;
            cursor: pointer;
        }
        </style>
        <body>
            %s
        </body>
    </html>"""
        html = html_string % (sigPvalue.to_html(
            classes='mystyle', escape=False))
        with open(file, 'w') as f:
            f.write(html)
        f.close()
        return html

    def saveModel(self, file='./learner.pkl'):
        """Save fitted learner

        Parameters
        ----------
        file : str, optional
            The path to save fitted learner, by default './learner.pkl'
        """
        with bz2.BZ2File(file + '.pbz2', 'w') as f:
            cPickle.dump(self, f)
        f.close()


class CircularLearner(BaseLearner, Circular):
    """Circular fragment leanrner

    Parameters
    ----------
    BaseLearner : smash.BaseLearner
        The base learner
    Circular : smash.fragments.Circular
        The class to obtain circular fragment matrix
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
        Circular.__init__(self, minRadius=minRadius, maxRadius=maxRadius,
                          nBits=nBits, folded=folded,
                          maxFragment=maxFragment, nJobs=nJobs)

    def GetMatrix(self, mols, **kwgrs):
        """Rewrite method

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            compounds, which have aspecific endpoint label, used to obtain significant fragments.

        Returns
        -------
        maxtrix : pandas.core.frame.DataFrame
            The calculated fragment matrix
        """
        matrix = self.GetCircularMatrix(mols, **kwgrs)
        return matrix

    def ShowFragment(self, mol, fragmentIndex):

        smarts, svg = self.ShowCircularFragment(mol, fragmentIndex)
        return (smarts, svg)

class PathLearner(BaseLearner, Path):
    """Path-Based fragment leanrner

    Parameters
    ----------
    BaseLearner : smash.BaseLearner
        The base learner
    Path : smash.fragments.Path
        The class to obtain path-based fragment matrix
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
        Path.__init__(self, minPath=minPath, maxPath=maxPath,
                      nBits=nBits, folded=folded, nJobs=nJobs, maxFragment=maxFragment)

    def GetMatrix(self, mols, **kwgrs):
        """Rewrite method

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            compounds, which have aspecific endpoint label, used to obtain significant fragments.

        Returns
        -------
        maxtrix : pandas.core.frame.DataFrame
            The calculated fragment matrix
        """
        matrix = self.GetPathMatrix(mols, **kwgrs)
        return matrix

    def ShowFragment(self, mol, fragmentIndex):
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
        smarts, svg = self.ShowPathFragment(mol, fragmentIndex)
        return smarts, svg


class FunctionGroupLearner(BaseLearner, FunctionGroup):
    """Function-Group fragment leanrner

    Parameters
    ----------
    BaseLearner : smash.BaseLearner
        The base learner
    FunctionGroup : smash.fragments.FunctionGroup
        The class to obtain Function-Group fragment matrix
    """

    def __init__(self, nJobs=1):
        """Initialization

        Parameters
        ----------
        nJobs : int, optional
            The number of CPUs to use to do the computation, by default 1
        """
        FunctionGroup.__init__(self, nJobs)

    def GetMatrix(self, mols, **kwgrs):
        """Rewrite method

        Parameters
        ----------
        mols : Iterable object, and each element is a rdkit.Chem.rdchem.Mol object.
            compounds, which have aspecific endpoint label, used to obtain significant fragments.

        Returns
        -------
        maxtrix : pandas.core.frame.DataFrame
            The calculated fragment matrix
        """
        matrix = self.GetFunctionGroupsMatrix(mols, **kwgrs)
        return matrix


if '__main__' == __name__:
    from rdkit import Chem
    from itertools import compress
    import time

    data = pd.read_csv(r'..\tests\carc\carc.txt', sep='\t')
    # data = pd.read_csv(r'..\tests\agg\agg.csv', sep=',')
#    data = data.sample(n=10000, replace=False)

    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    bo = mols.notna().values

    mols = list(compress(mols, bo))
    y_true = data.Label.values[bo]

    # start = time.clock()
    # circular = CircularLearner(minRadius=1, maxRadius=6,
    #                            maxFragment=True, nJobs=20)

    # sigPvalue, sigMatrix = circular.fit(mols, y_true)
    # sigPvalue.to_html('././092220_1.html', escape=False)
    # print(type(sigPvalue))
    # circular.savePvalue(sigPvalue, './092220.html')
    # end = time.clock()
    # print(end-start)

    start = time.clock()
    pa = PathLearner(
        minPath=1, maxPath=7,
        maxFragment=True, nJobs=20
    )

    sigPvalue, sigMatrix = pa.fit(mols, y_true)
    print(type(sigPvalue))
    pa.savePvalue(sigPvalue, './092220.html')
    end = time.clock()
    print(end-start)