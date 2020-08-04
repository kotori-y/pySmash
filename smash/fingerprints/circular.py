# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:13:34 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from rdkit import Chem
import pandas as pd
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import _getMorganEnv


__all__ = ['GetFoldedCircularFragment',
           'GetUnfoldedCircularFragment',
           'GetCircularFragment']


def _DisposeCircularBitInfo(mol, bitInfo, minRadius=3, maxFragment=True, svg=False):
    """dispose the bitinfo retrived from
    GetFoldedCircularFragment() or GetUnfoldedCircularFragment()

    Parameters
    ----------
    bitInfo : dict
        the key of dict is the number of bit (count from 0),
        and value is a tuple, the first one is
        the ceter atom, the second is the radius.
    maxFragment : bool, optional
        whether only return the max fragment of each atom, by default False

    Returns
    -------
    dict
        disposed bitinfo
    """
    idxAll = list(bitInfo.keys())

    if maxFragment:
        bitInfo = [[atom[0], atom[1], idx]
                   for idx, atomInfo in bitInfo.items() for atom in atomInfo]
        bitInfo = pd.DataFrame(bitInfo).groupby(0).max().drop_duplicates(2)
        idxFrag = bitInfo[2].tolist()
        bitInfo = tuple(zip(bitInfo.index, bitInfo[1]))

    else:
        idxFrag = idxAll
        bitInfo = tuple(map(lambda x: x[0], bitInfo.values()))

    fragments = {}
    for info, idx in zip(bitInfo, idxFrag):
        a, r = info

        if r >= minRadius:
            smi, svgImg = DrawMorganEnv(mol, a, r)
            fragments[idx] = (smi, svgImg.replace('\n','')) if svg else smi
        else:
            pass

    return (idxAll, fragments)


def DrawMorganEnv(mol, atomId, radius, molSize=(150, 150), baseRad=0.3, useSVG=True,
                  aromaticColor=(0.9, 0.9, 0.2), ringColor=(0.8, 0.8, 0.8),
                  centerColor=(0.6, 0.6, 0.9), extraColor=(0.9, 0.9, 0.9), drawOptions=None,
                  **kwargs):
    menv = _getMorganEnv(mol, atomId, radius, baseRad, aromaticColor, ringColor, centerColor,
                         extraColor, **kwargs)

    submol = menv.submol
    subMol = Chem.MolToSmiles(submol, canonical=True, isomericSmiles=False)
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


def GetFoldedCircularFragment(mol,
                              minRadius=3, maxRadius=6, 
                              nBits=1024, maxFragment=True, 
                              svg=False):
    """Get folded circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    radius : int, optional
        the radius of circular fingerprints. The default is 2.
    nBits : int, optional
        the number of bit of morgan. The default is 1024.
    maxFragment : bool, optional
        whether only return the max fragment of each atom, by default False

    Returns
    -------
    dict
        the key of dict is the number of bit (count from 0), 
        and value is a tuple, the first one is 
        the ceter atom, the second is the radius.

    """
    bitInfo = {}
    fp = GetMorganFingerprintAsBitVect(mol,
                                       radius=maxRadius,
                                       nBits=nBits,
                                       bitInfo=bitInfo)

    fragments = _DisposeCircularBitInfo(mol, bitInfo, minRadius, maxFragment, svg)
    return fragments


def GetUnfoldedCircularFragment(mol, minRadius=3, maxRadius=6, maxFragment=True, svg=False):
    """Get unfolded circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    radius : int, optional
        the radius of circular fingerprints. The default is 2.
    maxFragment : bool, optional
        whether only return the max fragment of each atom, by default False

    Returns
    -------
    bi : dict
        the key of dict is the number of bit (count from 0), 
        and value is a tuple of tuple, in each tuple element the first one is 
        the ceter atom, the second is the radius.

    """
    bitInfo = {}
    fp = GetMorganFingerprint(mol,
                              radius=maxRadius,
                              bitInfo=bitInfo)

    fragments = _DisposeCircularBitInfo(mol, bitInfo, minRadius, maxFragment, svg)
    return fragments


def GetCircularFragment(mol,
                        maxRadius=6, minRadius=3,
                        nBits=1024, folded=False,
                        maxFragment=True, svg=False):
    """Get circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    maxRadius : int, optional
        the maximum radius of circular fragment, by default 2
    minRadius : int, optional
        the minimum radius of circular fragment, by default 1
    nBits : int, optional
        the deminision of fragment folded, by default 1024,
        this parameter would be ignored if folded set as False
    folded : bool, optional
        whether hash the fragment, by default True
    maxFragment : bool, optional
        whether only return the max fragment of each atom, by default False

    Returns
    -------
    substrcutures : list
        circular fragment with specific radius
    """
    if folded:
        fragments = GetFoldedCircularFragment(mol,
                                              minRadius=minRadius, maxRadius=maxRadius,
                                              nBits=nBits, maxFragment=maxFragment,
                                              svg=svg)
    else:
        fragments = GetUnfoldedCircularFragment(mol,
                                                minRadius=minRadius, maxRadius=maxRadius,
                                                maxFragment=maxFragment, svg=svg)

    return fragments


if '__main__' == __name__:
    mol = Chem.MolFromSmiles('C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1')

    fragments = GetCircularFragment(
        mol, maxRadius=3, minRadius=1, maxFragment=False, svg=True)
    print(fragments[1])



