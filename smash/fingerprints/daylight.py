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



def CalculateDaylight(mol, minPath=1, maxPath=7, nBits=2048):
    """
    Calculate fingerprint and return the info of each bit

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    minPath : int, optional
        minimum number of bonds to include in the subgraphs. The default is 1.
    maxPath : int, optional
        maximum number of bonds to include in the subgraphs. The default is 1.
    nBits : int, optional
        the number of bit of morgan. The default is 2048.

    Returns
    -------
    bi : dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture
    fp : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        fingertprint.

    """
    bitInfo = {}
    fp = RDKFingerprint(mol, 
                        minPath=minPath, 
                        maxPath=maxPath,
                        fpSize=nBits,
                        bitInfo=bitInfo)
    return bitInfo

    
def CalculateSparseDaylight(mol, minPath=1, maxPath=7):
    """
    Calculate fingerprint and return the info of each bit

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    minPath : int, optional
        minimum number of bonds to include in the subgraphs. The default is 1.
    maxPath : int, optional
        maximum number of bonds to include in the subgraphs. The default is 1.

    Returns
    -------
    bi : dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture
    fp : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        fingertprint.

    """
    bitInfo = {}
    fp = UnfoldedRDKFingerprintCountBased(mol, 
                                          minPath=minPath, 
                                          maxPath=maxPath, 
                                          bitInfo=bitInfo)
    return bitInfo


def getBeginEndAtom(mol):
    head = {}
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        head[idx] = [begin, end]
    return head

def SmashMolWithDaylight(mol, 
                         minPath=1, 
                         maxPath=7, 
                         nBits=2048,
                         sparse=False):
    """
    

    Parameters
    ----------
    mol : TYPE
        DESCRIPTION.
    minPath : TYPE, optional
        DESCRIPTION. The default is 1.
    maxPath : TYPE, optional
        DESCRIPTION. The default is 7.
    nBits : TYPE, optional
        DESCRIPTION. The default is 2048.

    Returns
    -------
    substrcuture : dict
        the substruct of mol and frequency.

    """
    if sparse:
        bitInfo = CalculateSparseDaylight(mol,
                                          minPath=minPath,
                                          maxPath=maxPath,
                                          )
    else:
        bitInfo = CalculateDaylight(mol, 
                                    minPath=minPath,
                                    maxPath=maxPath,
                                    nBits=nBits,
                                    )
    
    substrcutures = []
    head = getBeginEndAtom(mol)
    
    for info in bitInfo.values():
        atomsToUse = [head[bond] for bond in info[0]]
        atomsToUse = set(sum(atomsToUse, []))
        smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=info[0])
        substrcutures.append(smi)
        
    substrcutures = list(set(substrcutures))
    return substrcutures


if '__main__' == __name__:
    from rdkit import Chem
    
    mol = Chem.MolFromSmiles('CNCC(O)c1ccc(O)c(O)c1')
    substrcutures_1 = SmashMolWithDaylight(mol, sparse=False)
    substrcutures_2 = SmashMolWithDaylight(mol, sparse=True)
    
    
    