import numpy as np

class Atom(dict):
    def __init__(self,):
        # self["style"] = {}
        pass

    def dist(self, other):
        return ((self["x"] - other["x"])**2 + (self["y"] - other["y"])**2 + (self["z"] - other["z"])**2)**0.5
    
    def set_style(self, style:dict):
        self.add_style(style)
        # self["style"] = style

    def add_style(self, style:dict):
        self["style"] = {**self["style"], **style} if "style" in self else style

    def get_style(self):
        return self["style"] if "style" in self else {}
    
    def remove_style(self):
        self.pop("style", None)



class Molecule(list):
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def __getitem__(self, idx):
        return [x for x in self.atoms if int(x['idx']) == int(idx)][0]
    
    def __iter__(self):
        for elem in self.atoms:
            yield elem
    
    def __len__(self):
        return len(self.atoms)
    
    def select_distance(self, atom: Atom, distance):
        return [a for a in self.atoms if atom.dist(a) <= distance]
    
    def get_mean(self):
        mean_x = sum([a["x"] for a in self.atoms])/len(self.atoms)
        mean_y = sum([a["y"] for a in self.atoms])/len(self.atoms)
        mean_z = sum([a["z"] for a in self.atoms])/len(self.atoms)
        return {"x": mean_x, "y": mean_y, "z": mean_z}

    def set_style(self, style:dict):
        for at in self.atoms:
            at.set_style(style)

    def get_coords(self):
        return np.asarray([[a["x"], a["y"], a["z"]] for a in self.atoms])
    
    def get_elements(self):
        return np.asarray([a["element"] for a in self.atoms])