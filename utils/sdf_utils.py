from .structures import Atom, Molecule

class SDFAtom(Atom):
    def __init__(self, i, line):
        super().__init__()

        self["idx"] = i
        self["x"] = float(line[0:10].strip())
        self["y"] = float(line[10:20].strip())
        self["z"] = float(line[20:30].strip())
        self["element"] = line[30:32].strip()
        self["remainder"] = line[32:]

    def __str__(self):
        line = list(" " * 69)
        line[0:10] = f'{self["x"]:.4f}'.rjust(10)
        line[10:20] = f'{self["y"]:.4f}'.rjust(10)
        line[20:30] = f'{self["z"]:.4f}'.rjust(10)
        line[31] = self["element"]
        line[32:] = self["remainder"]
        return "".join(line) + "\n"

class SDFBond(dict):
    def __init__(self, idx, line):
        self["from"] = line[0:3].strip()
        self["to"] = line[4:6].strip()
        self["type"] = line[6:9].strip()
        self["remainder"] = line[9:]
        self["idx"] = idx

    def __str__(self):
        line = list(" " * 11)
        line[0:3] = self["from"].rjust(3)
        line[3:5] = self["to"].rjust(3)
        line[6:9] = self["type"].rjust(3)
        line[9:12] = self["remainder"].rjust(3)
        return "".join(line) + "\n"
    
    def dist(self, other):
        return ((self["x"] - other["x"])**2 + (self["y"] - other["y"])**2 + (self["z"] - other["z"])**2)**0.5
    
    
class SDFMolecule(Molecule):
    def __init__(self, file):
        super().__init__()

        lines = file.split('\n')
        self.header = '\n'.join(lines[:2]) + '\n'
        self.count_line = lines[3]
        self.number_atoms = int(self.count_line.split()[0])
        self.number_bonds = int(self.count_line.split()[0])
        self.atoms = [SDFAtom(i, line) for i, line in enumerate(lines[4:4+self.number_atoms])]

        self.remainder = lines[4+self.number_atoms:]
        # self.bonds = [SDFBond(i, line) for i, line in enumerate(lines[4+self.number_atoms:4+self.number_atoms+self.number_bonds-1])]
        # self.remainder = lines[4+self.number_atoms+self.number_bonds:]

    def __str__(self):
        outstr = self.header + '\n'
        outstr += self.count_line + '\n'
        for at in self.atoms:
            outstr += str(at)
        for bond in self.bonds:
            outstr += str(bond)
        outstr += "\n".join(self.remainder)
        return outstr
    

    