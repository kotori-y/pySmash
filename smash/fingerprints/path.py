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

__all__ = ['GetFoldedPathFragment',
           'GetUnfoldedPathFragment',
           'GetPathFragment']


def _DisposePathFragments(fragments, maxFragment=True):

    def func(frags):              
        sets=[set(e) for e in frags]
        us=[]
        for e in sets:
            if any(e < s for s in sets):
                continue
            else:
                us.append(str(e))
        return us
    
    keys = list(fragments.keys())
    if maxFragment:
        dic = {str(set(x)):k for k,v in fragments.items() for x in v}
        fragments = [i for j in fragments.values() for i in j]
        fragments = func(fragments)
        fragments = {dic[k]:k for k in fragments}
    else:
        pass
    
    return (keys, fragments)





def GetFoldedPathFragment(mol, minPath=1, maxPath=7, nBits=2048, maxFragment=True):
    """Calculate folded path fragment.

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
    dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture

    """
    bitInfo = {}
    fp = RDKFingerprint(mol,
                        minPath=minPath,
                        maxPath=maxPath,
                        fpSize=nBits,
                        bitInfo=bitInfo)
    
    fragments = _DisposePathFragments(bitInfo, maxFragment=maxFragment)                  
    return fragments


def GetUnfoldedPathFragment(mol, minPath=1, maxPath=7, maxFragment=True):
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
    dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture

    """
    bitInfo = {}
    fp = UnfoldedRDKFingerprintCountBased(mol,
                                          minPath=minPath,
                                          maxPath=maxPath,
                                          bitInfo=bitInfo)

    fragments = _DisposePathFragments(bitInfo, maxFragment=maxFragment)                  
    return fragments


def getBeginEndAtom(mol):
    head = {}
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        head[idx] = [begin, end]
    return head


def GetPathFragment(mol,
                    minPath=1,
                    maxPath=7,
                    nBits=2048,
                    folded=False,
                    maxFragment=True):
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
    dict
        the substruct of mol and frequency.

    """
    if folded:
        bitInfo = GetFoldedPathFragment(mol,
                                        minPath=minPath,
                                        maxPath=maxPath,
                                        nBits=nBits,
                                        maxFragment=maxFragment)
    else:
        bitInfo = GetUnfoldedPathFragment(mol,
                                          minPath=minPath,
                                          maxPath=maxPath,
                                          maxFragment=maxFragment)

    # substrcutures = []
    # head = getBeginEndAtom(mol)

    # for info in bitInfo.values():
    #     atomsToUse = [head[bond] for bond in info[0]]
    #     atomsToUse = set(sum(atomsToUse, []))
    #     smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=info[0])
    #     substrcutures.append(smi)

    # substrcutures = list(set(substrcutures))
    # return substrcutures
    return bitInfo

if '__main__' == __name__:
    
    mol = Chem.MolFromSmiles('CN(C1=CC=C(C=C1)N=N/C2=CC=CC=C2)C')
    frag = GetPathFragment(mol, maxFragment=True)
    print(frag)
