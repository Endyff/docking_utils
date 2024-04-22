from .structures import Atom, Molecule

class PDBAtom(Atom):
    def __init__(self, line):
        super().__init__()

        self["type"] = line[0:6].strip()
        self["idx"] = line[6:11].strip()
        self["name"] = line[12:16].strip()
        self["resn"] = line[17:20].strip()
        self["resid"] = int(int(line[22:26]))
        self["x"] = float(line[30:38])
        self["y"] = float(line[38:46])
        self["z"] = float(line[46:54])
        self["element"] = line[76:78].strip()

    def __str__(self):
        line = list(" " * 80)
        line[0:6] = self["type"].ljust(6)
        line[6:11] = self["idx"].ljust(5)
        line[12:16] = self["name"].ljust(4)
        line[17:20] = self["resn"].ljust(3)
        line[22:26] = str(self["resid"]).ljust(4)
        line[30:38] = str(self["x"]).rjust(8)
        line[38:46] = str(self["y"]).rjust(8)
        line[46:54] = str(self["z"]).rjust(8)
        line[76:78] = self["element"].rjust(2)
        return "".join(line) + "\n"
    

class PDBMolecule(Molecule):
    def __init__(self, file):
        super().__init__()
        self.atoms = [PDBAtom(line) for line in file.split('\n') if "ATOM" == line[:4]]
    def __str__(self):
        outstr = ""
        for at in self.atoms:
            outstr += str(at)

        return outstr
    
    def select_residue(self, atom: PDBAtom):
        return [a for a in self.atoms if a["resn"] == atom["resn"] and a["resid"] == atom["resid"]]
    
    def set_residue_style(self, atom:PDBAtom, style:dict):
        for at in self.select_residue(atom):
            at.add_style(style)