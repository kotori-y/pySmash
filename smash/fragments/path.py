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
           'GetPathFragment']


def _DisposePathFragments(mol, bitInfo, maxFragment=True, svg=False):
    """Dispose the bitinfo retrived from GetFoldedPathFragment() or GetUnfoldedPathFragment()
    *internal only*
    """
    def func(frags):
        sets = [set(e) for e in frags]
        us = []
        for e in sets:
            if any(e < s for s in sets):
                continue
            else:
                us.append(str(e))
        return us

    fragments = {}
    keys = list(bitInfo.keys())
    if maxFragment:
        dic = {str(set(x)): k for k, v in bitInfo.items() for x in v}
        bitInfo = [i for j in bitInfo.values() for i in j]
        bitInfo = func(bitInfo)
        for k in bitInfo:
            bondPath = list(eval(k))
            smi, svgImg = DrawRDKitEnv(mol, bondPath)
            fragments[dic[k]] = (smi, svgImg.replace('\n','')) if svg else (smi, )
    else:
        pass

    return (keys, fragments)


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
    return subMol, drawer.GetDrawingText()


def GetFoldedPathFragment(mol,
                          minPath=1, maxPath=7,
                          nBits=1024, maxFragment=True,
                          svg=False):
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
    svg : bool, optional
        Whether output with a svg image, by default False

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
        mol, bitInfo, maxFragment=maxFragment, svg=svg)
    return fragments


def GetUnfoldedPathFragment(mol, minPath=1, maxPath=7,
                            maxFragment=True, svg=False):
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
    svg : bool, optional
        Whether output with a svg image, by default False

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

    fragments = _DisposePathFragments(mol, bitInfo, maxFragment=maxFragment, svg=svg)
    return fragments


def GetPathFragment(mol,
                    minPath=1, maxPath=7,
                    nBits=1024, folded=False,
                    maxFragment=True, svg=False):
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
    svg : bool, optional
        Whether output with a svg image, by default False

    Returns
    -------
    fragments : tuple
        The first element is the ID of all fragments generated,
        and the second one is a dict whose key is the ID of output fragments,
        value is corresponding SMARTS and svg string (is svg set as True)
    """
    if folded:
        fragments = GetFoldedPathFragment(mol,
                                        minPath=minPath, maxPath=maxPath,
                                        nBits=nBits, maxFragment=maxFragment,
                                        svg=svg)
    else:
        fragments = GetUnfoldedPathFragment(mol,
                                          minPath=minPath, maxPath=maxPath,
                                          maxFragment=maxFragment, svg=svg)

    # substrcutures = []
    # head = getBeginEndAtom(mol)

    # for info in bitInfo.values():
    #     atomsToUse = [head[bond] for bond in info[0]]
    #     atomsToUse = set(sum(atomsToUse, []))
    #     smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=info[0])
    #     substrcutures.append(smi)

    # substrcutures = list(set(substrcutures))
    # return substrcutures
    return fragments


if '__main__' == __name__:

    mol = Chem.MolFromSmiles('CN(C1=CC=C(C=C1)N=N/C2=CC=CC=C2)C')
    fragments = GetPathFragment(mol, maxFragment=True, svg=True)
    print(fragments)
