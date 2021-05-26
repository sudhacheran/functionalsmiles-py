from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
import re
from rdkit.Chem import rdMolHash
import util.rdkit_func as rkf
import sys

# static method to display the results
def generate_data_table(reactant, product):
    data_table = []
    out_title = ["", "Reactant", "Product"]
    data_table.append(out_title)
    fnsmi = FunctionSmiles()
    data_table.append(["Smiles", reactant, product])
    smi_tk1 = ["Regioisomers", fnsmi.rdkit_regioisomers(Chem.MolFromSmiles(reactant)), fnsmi.rdkit_regioisomers(Chem.MolFromSmiles(reactant))]
    data_table.append(smi_tk1)
    smi_tk2 = ["Func.SMILEs", fnsmi.get_functionalsmiles(reactant), fnsmi.get_functionalsmiles(product)]
    data_table.append(smi_tk2)
    return data_table


class FunctionSmiles:

    def __init__(self) -> None:
        super().__init__()

    def rdkit_regioisomers(self, mol):
        rg = rdMolHash.HashFunction.names.get("Regioisomer")
        regioisomers = rdMolHash.MolHash(mol, rg)
        return regioisomers

    def _fragment_on_bond(self,mol, atom1, atom2, posloc):
        bond = mol.GetBondBetweenAtoms(atom1, atom2)
        #print(atom1, atom2, bond)
        new_mol = mol
        if bond is not None:
            posp = int(str(posloc)+"0")
            posc = int(str(posloc)+"1")
            new_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], dummyLabels=[(posloc, posloc)])
        return new_mol

    def fragment(self, mol):
        new_mol = mol
        #print(mol)
        for i in range(mol.GetNumAtoms() - 1):
            atom1 = mol.GetAtomWithIdx(i)
            atom2 = mol.GetAtomWithIdx(i + 1)
            atom1_isaromatic = atom1.GetIsAromatic()
            atom2_isaromatic = atom2.GetIsAromatic()
            #print(atom1.GetSmarts(), atom2.GetSmarts())
            #print(atom1_isaromatic,  atom2_isaromatic)
            if (atom1_isaromatic and not atom2_isaromatic) or (atom2_isaromatic and not atom1_isaromatic):
                #print(atom1.GetSmarts(),atom2.GetSmarts())
                new_mol = self._fragment_on_bond(new_mol, atom1.GetIdx(), atom2.GetIdx(), self.posloc)
                self.posloc += 1
        return new_mol

    def tokenizeingroups(self,mol):
        prev_atom = ""
        for i in range(mol.GetNumAtoms()):
            atom1 = mol.GetAtomWithIdx(i)
            atom1_isaromatic = atom1.GetIsAromatic()
            prev_atom_isaromatic = False if prev_atom is "" else prev_atom.GetIsAromatic()
            print(prev_atom_isaromatic, atom1_isaromatic)
            if (atom1_isaromatic and not prev_atom_isaromatic):
                print(" ")
            print(atom1.GetSmarts())
            prev_atom = atom1

    def mol_with_atom_index(self,mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return mol

    def functionalsmileform(self,strsmile):
        strret = strsmile
        # strret= re.sub(r"\[[0-9]\*\]", "", strsmile)
        # print(strret)
        # strret2 =re.sub(r"\(([^)])+\)","",strret)
        # print(strret2)
        # strret =strret.replace("()","")
        # strret =strret.replace("*","")
        return strret

    def fragment2(self,mol):
        new_mol = mol
        tup1 = new_mol.GetSubstructMatches(Chem.MolFromSmarts('aA'))
        # tup2 = mol_tt.GetSubstructMatches(Chem.MolFromSmarts('A[OH]'))
        tup3 = mol.GetSubstructMatches(Chem.MolFromSmarts('[A;!R][A;R]'))
        tupgrp = tup1 + tup3  # + tup2

        for data in tupgrp:
            atom1 = mol.GetAtomWithIdx(data[0])
            atom2 = mol.GetAtomWithIdx(data[1])
            new_mol = self._fragment_on_bond(new_mol, atom1.GetIdx(), atom2.GetIdx(), self.posloc)
            self.posloc += 1
        return new_mol

    def groupdata(self, mol):
        if type(mol) is str:
            listmol2 = mol.split(".")
        else:
            listmol2 = Chem.MolToSmiles(mol).split(".")
        fsmi = " ".join(listmol2)
        return fsmi

    def get_functionalsmiles(self, smi):
        self.posloc = 1
        self.smi = smi
        if self.smi.__contains__("."):
            self.smi_list = self.smi.split(".")
            self.cp_mol = [Chem.MolFromSmiles(smistr) for smistr in self.smi_list]
        else:
            self.cp_mol = Chem.MolFromSmiles(smi)

        if type(self.cp_mol) is list:
           moldata = [self.fragment2(mol1) for mol1 in self.cp_mol]
           retsmi = " . ".join(self.groupdata(mol2) for mol2 in moldata)
        else:
           moldata = self.fragment2(self.cp_mol)
           retsmi = self.groupdata(moldata)
        return retsmi

    def findpos(self,smi):  # This includes {} and () formatting
        pattern = "\[([0-9]*?\*)\]"
        regex = re.compile(pattern)
        tokens = ['[' + token + ']' for token in regex.findall(smi)]
        return tokens

    def regroupfunctioalsmiles(self, strSmi):
        compounds = strSmi.split(".")
        compounds = [smi.strip() for smi in compounds]
        #print("Fragmented", compounds)
        posremoved = []
        compound2 = []
        fragmol = ""
        for molecule in compounds:
            if "*" in molecule:
                fragmol = ".".join(molecule.split(" "))
                pos = self.findpos(fragmol)
                for mks in pos:
                    if mks not in posremoved:
                        edmol = Chem.MolFromSmarts(fragmol) if type(fragmol) is str else fragmol
                        # print(Chem.MolToSmiles(edmol))
                        edcombo = Chem.EditableMol(edmol)
                        # print(mks)
                        posremoved.append(mks)
                        ss = edmol.GetSubstructMatches(Chem.MolFromSmarts("*" + mks))
                        # print(ss)
                        edcombo.AddBond(ss[0][0], ss[1][0], order=Chem.rdchem.BondType.SINGLE)
                        fragmol = edcombo.GetMol()
                        fragmol = AllChem.DeleteSubstructs(fragmol, Chem.MolFromSmarts(mks))
                try:
                    Chem.SanitizeMol(fragmol)
                except Chem.rdchem.KekulizeException:
                   print(sys.exc_info()[0])
                   pass
                compound2.append(fragmol)
            else:
                compound2.append(Chem.MolFromSmiles(molecule))
        finalcmp = ".".join(Chem.MolToSmiles(cmp) for cmp in compound2)
        return finalcmp


if __name__ == "__main__":

    reactantsmile = "C C O . C C O C ( = O ) C N 1 C C C ( N 2 C C N ( C ( = O ) C ( C c 3 c c ( C ) c ( O C c 4 c c c c c 4 ) c ( C ) c 3 ) O C ( = O ) N 3 C C C ( N 4 C C c 5 c c c c c 5 N C 4 = O ) C C 3 ) C C 2 ) C C 1 . [H] [H]"
    productsmile = "COC(=O)c1ccccc1-c1ccc(COc2cc(C)nc3cccnc23)cc1"
    #print(generate_data_table(reactantsmile,productsmile))
    fnsmi = FunctionSmiles()
    reactantsmile = rkf.smi_detokenize(reactantsmile)
    print(reactantsmile)
    #print(Chem.CanonSmiles(reactantsmile))
    funcSmi=fnsmi.get_functionalsmiles(reactantsmile)
    print(funcSmi)
    print(fnsmi.regroupfunctioalsmiles(funcSmi))
    #print(fnsmi.get_functionalsmiles(productsmile))