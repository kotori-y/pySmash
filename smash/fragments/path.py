# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:13:29 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from rdkit import Chem
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.rdmolops import UnfoldedRDKFingerprintCountBased
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import _getRDKitEnv

__all__ = ['GetFoldedPathFragment',
           'GetUnfoldedPathFragment',
           'GetPathFragment',
           'DrawRDKitEnv']


def _DisposePathFragments(mol, bitInfo, maxFragment=True):
    """Dispose the bitinfo retrived from GetFoldedPathFragment() or GetUnfoldedPathFragment()
    *internal only*
    """
    allFragments = list(bitInfo.keys())
    st = []
    maxFragments = []
    
    for idx in sorted(bitInfo, key=lambda idx: len(bitInfo[idx][0]), reverse=True):
        smi, _ = DrawRDKitEnv(mol, bitInfo[idx][0])
        newPatt = Chem.MolFromSmarts(smi)
        newPatt.UpdatePropertyCache(strict=False)
                
        for patt in st:
            if not maxFragment or patt.HasSubstructMatch(newPatt):
                break
        else:     
            st.append(newPatt)
            maxFragments.append(idx)

    return (allFragments, maxFragments)


def DrawRDKitEnv(mol, bondPath, molSize=(150, 150), baseRad=0.3, useSVG=True,
                 aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9, 0.9), nonAromaticColor=None,
                 drawOptions=None, **kwargs):
    """Get SMARTS and SVG image from given bonds
    *internal only*
    """
    menv = _getRDKitEnv(mol, bondPath, baseRad, aromaticColor,
                        extraColor, nonAromaticColor, **kwargs)
    submol = menv.submol
    subMol = Chem.MolToSmiles(submol, isomericSmiles=False)

    # Drawing
    if useSVG:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])

    if drawOptions is None:
        drawopt = drawer.drawOptions()
        drawopt.continuousHighlight = False
    else:
        drawOptions.continuousHighlight = False
        drawer.SetDrawOptions(drawOptions)

    drawer.DrawMolecule(menv.submol, highlightAtoms=menv.highlightAtoms,
                        highlightAtomColors=menv.atomColors, highlightBonds=menv.highlightBonds,
                        highlightBondColors=menv.bondColors, highlightAtomRadii=menv.highlightRadii,
                        **kwargs)
    drawer.FinishDrawing()
    return subMol, drawer.GetDrawingText().replace('\n','')


def GetFoldedPathFragment(mol,
                          minPath=1, maxPath=7,
                          nBits=1024, maxFragment=True,
                          disposed=True):
    """Calculate folded path fragment.

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
    minPath : int, optional
        The probable minimum length of path-based fragment, by default 1
    maxPath : int, optional
        The probable maximum length of path-based fragment, by default 7
    nBits : int, optional
        the number of bit of morgan, by default 1014
    maxFragment : bool, optional
        Whether only return the maximum fragment of a given start atom, by default True
    disposed : bool, optional
        Whether dispose the original bitinfo, by default True

    Returns
    -------
    fragments : tuple
        The first element is the ID of all fragments generated,
        and the second one is a dict whose key is the ID of output fragments,
        value is corresponding SMARTS and svg string (is svg set as True)
    """
    bitInfo = {}
    fp = RDKFingerprint(mol,
                        minPath=minPath,
                        maxPath=maxPath,
                        fpSize=nBits,
                        bitInfo=bitInfo)

    fragments = _DisposePathFragments(
            mol, bitInfo, maxFragment=maxFragment
        ) if disposed else bitInfo
    return fragments


def GetUnfoldedPathFragment(mol, minPath=1, maxPath=7,
                            maxFragment=True, disposed=True):
    """Calculate unfolded path fragment.

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
    minPath : int, optional
        The probable minimum length of path-based fragment, by default 1
    maxPath : int, optional
        The probable maximum length of path-based fragment, by default 7
    maxFragment : bool, optional
        Whether only return the maximum fragment of a given start atom, by default True
    disposed : bool, optional
        Whether dispose the original bitinfo, by default True

    Returns
    -------
    fragments : tuple
        The first element is the ID of all fragments generated,
        and the second one is a dict whose key is the ID of output fragments,
        value is corresponding SMARTS and svg string (is svg set as True)
    """
    bitInfo = {}
    fp = UnfoldedRDKFingerprintCountBased(mol,
                                          minPath=minPath,
                                          maxPath=maxPath,
                                          bitInfo=bitInfo)

    fragments = _DisposePathFragments(
            mol, bitInfo, maxFragment=maxFragment
        ) if disposed else bitInfo
    return fragments


def GetPathFragment(mol,
                    minPath=1, maxPath=7,
                    nBits=1024, folded=False,
                    maxFragment=True, disposed=True):
    """Calculate path fragment.

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
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
    disposed : bool, optional
        Whether dispose the original bitinfo, by default True

    Returns
    -------
    fragments : tuple
        The first element is the ID of all fragments generated,
        and the second one is a dict whose key is the ID of output fragments,
        value is corresponding SMARTS and svg string (is svg set as True)
    """
    if folded:
        fragments = GetFoldedPathFragment(
            mol, minPath=minPath, maxPath=maxPath,
            nBits=nBits, maxFragment=maxFragment,
            disposed=disposed
        )
    else:
        fragments = GetUnfoldedPathFragment(
            mol, minPath=minPath, maxPath=maxPath,
            maxFragment=maxFragment, disposed=disposed
        )

    return fragments


if '__main__' == __name__:

    mol = Chem.MolFromSmiles('CN(C1=CC=C(C=C1)N=N/C2=CC=CC=C2)C')
    fragments = GetPathFragment(mol, maxFragment=True, disposed=True)
    print(fragments)
