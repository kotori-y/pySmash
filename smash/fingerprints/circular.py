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
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint

__all__ = ['GetFoldedCircularFragment',
           'GetUnfoldedCircularFragment',
           'GetCircularFragment']


def GetFoldedCircularFragment(mol, radius=2, nBits=1024):
    """Get folded circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    radius : int, optional
        the radius of circular fingerprints. The default is 2.
    nBits : int, optional
        the number of bit of morgan. The default is 1024.

    Returns
    -------
    bi : dict
        the key of dict is the number of bit (count from 0), 
        and value is a tuple of tuple, in each tuple element the first one is 
        the ceter atom, the second is the radius.

    """
    biInfo = {}
    fp = GetMorganFingerprintAsBitVect(mol,
                                       radius=radius,
                                       nBits=nBits,
                                       bitInfo=biInfo,
                                       )
    return biInfo


def GetUnfoldedCircularFragment(mol, radius=2):
    """Get unfolded circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    radius : int, optional
        the radius of circular fingerprints. The default is 2.

    Returns
    -------
    bi : dict
        the key of dict is the number of bit (count from 0), 
        and value is a tuple of tuple, in each tuple element the first one is 
        the ceter atom, the second is the radius.

    """
    bitInfo = {}
    fp = GetMorganFingerprint(mol,
                              radius=radius,
                              bitInfo=bitInfo,
                              )
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
    #finalsmi = copy.deepcopy(smi)
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
                        radius=2,
                        minRadius=1,
                        nBits=1024,
                        folded=True,
                        ):
    """Get circular fragment under specific radius

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    radius : int, optional
        the maximum radius of circular fragment, by default 2
    minRadius : int, optional
        the minimum radius of circular fragment, by default 1
    nBits : int, optional
        the deminision of fragment folded, by default 1024,
        this parameter would be ignored if folded set as False
    folded : bool, optional
        whether hash the fragment, by default True

    Returns
    -------
    substrcutures : list
        circular fragment with specific radius
    """
    if folded:
        bitInfo = GetFoldedCircularFragment(mol,
                                            radius=radius,
                                            nBits=nBits)
    else:
        bitInfo = GetUnfoldedCircularFragment(mol,
                                              radius=radius)

    substrcutures = []

    for info in bitInfo.values():
        a, r = info[0]

        if r >= minRadius:
            smi2 = getSubstructSmi(mol, a, r)
            substrcutures.append(smi2)
        else:
            pass

    substrcutures = list(set(substrcutures))
    return substrcutures


if '__main__' == __name__:
    mol = Chem.MolFromSmiles('CNCC(O)c1ccc(O)c(O)c1')
    substrcutures = GetCircularFragment(
        mol, radius=2, minRadius=1, folded=True)
    print(len(substrcutures))
