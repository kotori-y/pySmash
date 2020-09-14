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


def _DisposeCircularBitInfo(mol, bitInfo, minRadius=1, maxFragment=True, svg=False):
    """Dispose the bitinfo retrived from GetFoldedCircularFragment() or GetUnfoldedCircularFragment()
    *internal only*
    """
    idxAll = list(bitInfo.keys())

    if maxFragment:
        bitInfo = [[atom[0], atom[1], idx]
                   for idx, atomInfo in bitInfo.items() for atom in atomInfo]
        bitInfo = pd.DataFrame(bitInfo, columns=['atom', 'radius', 'idx'])
        bitInfo = bitInfo.loc[bitInfo.groupby('atom').radius.idxmax()]
        idxFrag = bitInfo['idx'].to_list()
        bitInfo = tuple(zip(bitInfo.atom.to_list(), bitInfo.radius.to_list()))

    else:
        idxFrag = idxAll
        bitInfo = tuple(map(lambda x: x[0], bitInfo.values()))

    fragments = {}
    for info, idx in zip(bitInfo, idxFrag):
        a, r = info

        if r >= minRadius:
            smi, svgImg = DrawMorganEnv(mol, a, r)
            fragments[idx] = (smi, svgImg.replace(
                '\n', '')) if svg else (smi, )
        else:
            pass

    return (idxAll, fragments)


def DrawMorganEnv(mol, atomId, radius, molSize=(150, 150), baseRad=0.3, useSVG=True,
                  aromaticColor=(0.9, 0.9, 0.2), ringColor=(0.8, 0.8, 0.8),
                  centerColor=(0.6, 0.6, 0.9), extraColor=(0.9, 0.9, 0.9), drawOptions=None,
                  **kwargs):
    """Get SMARTS and SVG image from given ceten and radius
    *internal only*
    """
    menv = _getMorganEnv(mol, atomId, radius, baseRad, aromaticColor, ringColor, centerColor,
                         extraColor, **kwargs)

    submol = menv.submol
    subMol = Chem.MolToSmiles(submol, isomericSmiles=False, allBondsExplicit=True, allHsExplicit=True)
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


def GetFoldedCircularFragment(mol, minRadius=1, maxRadius=2,
                              nBits=1024, maxFragment=True,
                              svg=False):
    """Get folded circular fragment

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
    minRadius : int, optional
        The probable minimum radius of circular fragment, by default 1
    maxRadius : int, optional
        The probable maximum radius of circular fragment, by default 2
    nBits : int, optional, 
        the number of bit of morgan, by default 1014
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
    bitInfo = {}
    fp = GetMorganFingerprintAsBitVect(mol,
                                       radius=maxRadius,
                                       nBits=nBits,
                                       bitInfo=bitInfo)

    fragments = _DisposeCircularBitInfo(
        mol, bitInfo, minRadius, maxFragment, svg)
    return fragments


def GetUnfoldedCircularFragment(mol, minRadius=1, maxRadius=2,
                                maxFragment=True, svg=False):
    """Get unfolded circular fragment

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
    minRadius : int, optional
        The probable minimum radius of circular fragment, by default 1
    maxRadius : int, optional
        The probable maximum radius of circular fragment, by default 2
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
    bitInfo = {}
    fp = GetMorganFingerprint(mol,
                              radius=maxRadius,
                              bitInfo=bitInfo)

    fragments = _DisposeCircularBitInfo(
        mol, bitInfo, minRadius, maxFragment, svg)
    return fragments


def GetCircularFragment(mol, minRadius=1, maxRadius=2,
                        nBits=1024, folded=False,
                        maxFragment=True, svg=False):
    """Get circular fragment

    Parameters
    ----------
    mol : dkit.Chem.rdchem.Mol object
        Compound to be Calculated
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

    Returns
    -------
    fragments : tuple
        The first element is the ID of all fragments generated,
        and the second one is a dict whose key is the ID of output fragments,
        value is corresponding SMARTS and svg string (is svg set as True)
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
        mol, maxRadius=4, minRadius=1, maxFragment=False, svg=True)
    print(fragments)
