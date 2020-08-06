# -*- encoding: utf-8 -*-
'''
Created on 2020/07/03 21:46:19

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
'''


#  Original authors: Richard Hall and Guillaume Godin
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

#
#
# Richard hall 2017
# IFG main code
# Guillaume Godin 2017
# refine output function
# astex_ifg: identify functional groups a la Ertl, J. Cheminform (2017) 9:36


from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

__all__ = ['GetFunctionGroupFragment']


subMolDic = {}


def _DisposedFgFragment(mol, subMol, molSize=(150, 150)):
    if not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)

    subMolDic[subMol] = Chem.MolFromSmarts(subMol)
    atomToUsed = mol.GetSubstructMatch(subMolDic[subMol])
    amap = {}
    submol = Chem.PathToSubmol(mol, atomToUsed, atomMap=amap)

    Chem.FastFindRings(submol)
    conf = Chem.Conformer(submol.GetNumAtoms())
    confOri = mol.GetConformer(0)
    for i1, i2 in amap.items():
        conf.SetAtomPosition(i2, confOri.GetAtomPosition(i1))
    submol.AddConformer(conf)

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawopt = drawer.drawOptions()
    drawopt.continuousHighlight = False
    drawer.DrawMolecule(submol)
    drawer.FinishDrawing()

    return drawer.GetDrawingText().replace('\n', '')


def merge(mol, marked, aset):
    bset = set()
    for idx in aset:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            jdx = nbr.GetIdx()
            if jdx in marked:
                marked.remove(jdx)
                bset.add(jdx)
    if not bset:
        return
    merge(mol, marked, bset)
    aset.update(bset)


# atoms connected by non-aromatic double or triple bond to any heteroatom
# c=O should not match (see fig1, box 15).  I think using A instead of * should sort that out?
PATT_DOUBLE_TRIPLE = Chem.MolFromSmarts('A=,#[!#6]')
# atoms in non aromatic carbon-carbon double or triple bonds
PATT_CC_DOUBLE_TRIPLE = Chem.MolFromSmarts('C=,#C')
# acetal carbons, i.e. sp3 carbons connected to tow or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds
PATT_ACETAL = Chem.MolFromSmarts('[CX4](-[O,N,S])-[O,N,S]')
# all atoms in oxirane, aziridine and thiirane rings
PATT_OXIRANE_ETC = Chem.MolFromSmarts('[O,N,S]1CC1')

PATT_TUPLE = (PATT_DOUBLE_TRIPLE, PATT_CC_DOUBLE_TRIPLE,
              PATT_ACETAL, PATT_OXIRANE_ETC)


def GetFunctionGroupFragment(mol, svg=False):
    marked = set()
# mark all heteroatoms in a molecule, including halogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1):  # would we ever have hydrogen?
            marked.add(atom.GetIdx())

# mark the four specific types of carbon atom
    for patt in PATT_TUPLE:
        for path in mol.GetSubstructMatches(patt):
            for atomindex in path:
                marked.add(atomindex)

# merge all connected marked atoms to a single FG
    groups = []
    while marked:
        grp = set([marked.pop()])
        merge(mol, marked, grp)
        groups.append(grp)

# extract also connected unmarked carbon atoms
    ifgs = []
    for g in groups:
        uca = set()
        for atomidx in g:
            for n in mol.GetAtomWithIdx(atomidx).GetNeighbors():
                if n.GetAtomicNum() == 6:
                    uca.add(n.GetIdx())
        fg = Chem.MolFragmentToSmiles(mol, g.union(uca), canonical=True)
        if fg not in ifgs:
            if svg:
                ifgs.append((fg, _DisposedFgFragment(mol, fg)))
            else:
                ifgs.append((fg, ))
    return ifgs


if __name__ == "__main__":
    from rdkit import Chem
    from itertools import compress
    import pandas as pd

    data = pd.read_csv(r'tests\Canc\Canc.txt', sep='\t')
    data = data.sample(n=100)

    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    bo = mols.notna().values

    mols = list(compress(mols, bo))
    
    for mol in mols:
        fgs = GetFunctionGroupFragment(mol, svg=True)
