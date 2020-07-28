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

__all__ = ['GetFoldedCircularFragment',
           'GetUnfoldedCircularFragment',
           'GetCircularFragment']


def _DisposeCircularBitInfo(mol, bitInfo, minRadius=3, maxFragment=False):
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
            smi2 = getSubstructSmi(mol, a, r)
            fragments[idx] = smi2
        else:
            pass

    return (idxAll, fragments)


def GetFoldedCircularFragment(mol, minRadius=3, maxRadius=6, nBits=1024, maxFragment=False):
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

    bitInfo = _DisposeCircularBitInfo(mol, bitInfo, minRadius, maxFragment)
    return bitInfo


def GetUnfoldedCircularFragment(mol, minRadius=3, maxRadius=6, maxFragment=False):
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

    bitInfo = _DisposeCircularBitInfo(mol, bitInfo, minRadius, maxFragment)
    return bitInfo


def includeRingMembership(s, n):
    r = ';R]'
    d = "]"
    return r.join([d.join(s.split(d)[:n]), d.join(s.split(d)[n:])])


def includeDegree(s, n, d):
    r = ';D'+str(d)+']'
    d = "]"
    return r.join([d.join(s.split(d)[:n]), d.join(s.split(d)[n:])])


def writePropsToSmiles(mol, smi, order):
    # finalsmi = copy.deepcopy(smi)
    finalsmi = smi
    for i, a in enumerate(order):
        atom = mol.GetAtomWithIdx(a)
        if atom.IsInRing():
            finalsmi = includeRingMembership(finalsmi, i+1)
        finalsmi = includeDegree(finalsmi, i+1, atom.GetDegree())
    return finalsmi


def getSubstructSmi(mol, atomID, radius):
    if radius > 0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atomID)
        atomsToUse = []
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env = None

    smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=env,
                                   allHsExplicit=True, allBondsExplicit=True,
                                   rootedAtAtom=atomID)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    smi2 = writePropsToSmiles(mol, smi, order)
    return smi2


def GetCircularFragment(mol,
                        maxRadius=6,
                        minRadius=3,
                        nBits=1024,
                        folded=True,
                        maxFragment=False):
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
        bitInfo = GetFoldedCircularFragment(mol,
                                            minRadius=minRadius,
                                            maxRadius=maxRadius,
                                            nBits=nBits,
                                            maxFragment=maxFragment)
    else:
        bitInfo = GetUnfoldedCircularFragment(mol,
                                              minRadius=minRadius,
                                              maxRadius=maxRadius,
                                              maxFragment=maxFragment)

    # substrcutures = []

    # for info in bitInfo[1]:
    #     a, r = info

    #     if r >= minRadius:
    #         smi2 = getSubstructSmi(mol, a, r)
    #         substrcutures.append(smi2)
    #     else:
    #         pass

    # substrcutures = list(set(substrcutures))
    return bitInfo


if '__main__' == __name__:
    mol = Chem.MolFromSmiles('C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1')

    bitInfo = GetUnfoldedCircularFragment(
        mol, maxRadius=3, minRadius=1, maxFragment=False)
    print(bitInfo[1].keys())

    # fragments = GetCircularFragment(mol,
    #                                 folded=False, minRadius=3,
    #                                 maxRadius=6, maxFragment=True)
    # print(fragments)
