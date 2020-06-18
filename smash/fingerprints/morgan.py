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



def CalculateMorgan(mol, radius=2, nBits=1024):
    """
    Calculate morgan fingerprint and return the info of NonzeroElements

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
   

def CalculateSparseMorgan(mol, radius=2):
    """
    Calculate sparse morgan fingerprint and return the info of NonzeroElements

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
    r=';R]'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def includeDegree(s, n, d):
    r=';D'+str(d)+']'
    d="]"
    return r.join([d.join(s.split(d)[:n]),d.join(s.split(d)[n:])])
 
def writePropsToSmiles(mol, smi, order):
    #finalsmi = copy.deepcopy(smi)
    finalsmi = smi
    for i,a in enumerate(order):
        atom = mol.GetAtomWithIdx(a)
        if atom.IsInRing():
            finalsmi = includeRingMembership(finalsmi, i+1)
        finalsmi = includeDegree(finalsmi, i+1, atom.GetDegree())
    return finalsmi
 
def getSubstructSmi(mol, atomID, radius):
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env=None
    
    smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=env,
                                   allHsExplicit=True, allBondsExplicit=True, 
                                   rootedAtAtom=atomID)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    smi2 = writePropsToSmiles(mol,smi,order)
    return smi2

def SmashMolWithMorgan(mol, 
                       radius=2, 
                       minRadius=1,
                       nBits=1024, 
                       sparse=False,
                       ):
    """
    
    """
    if sparse:
        bitInfo = CalculateSparseMorgan(mol, 
                                        radius=radius)
    else:
        bitInfo = CalculateMorgan(mol, 
                                  radius=radius, 
                                  nBits=nBits)
    
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
    substrcutures = SmashMolWithMorgan(mol, radius=3, minRadius=2, sparse=True)
    print(substrcutures)