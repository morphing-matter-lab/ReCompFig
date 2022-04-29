"""
name = "Nylon"
poisson_ratio = 0.42
elastic_mod = 1800 # MPA
shear_mod = elastic_mod / (2 * (1 + poisson_ratio))
max_strain = 0.15 # %
density = 1.01 # g/cm3

name = "TPU"
poisson_ratio = 0.3897
elastic_mod = 60 # MPA
shear_mod = elastic_mod / (2 * (1 + poisson_ratio))
max_strain = 0.9 # %
density = 1.11 # g/cm3
"""

class Material(object):
    def __init__(self, id, name):
        self.id = id
        self.name = name

    @staticmethod
    def create(info):

        objType = info[0]
        for elemType in Material.__subclasses__():
            if objType == elemType.name:
                modeled = elemType(info)

        return modeled

class IsoMaterial(Material):
    name = "Isotropic"

    def __init__(self, info):
        super().__init__(info[1], info[2])

        self.E = info[3] # elastic modulus
        self.G = info[4] # shear modulus
        self.D = info[5] # density
        self.max_strain = info[6] # max strain

    def __repr__(self):

        msg = "Isotropic material(ID:%d, name:%s)" % (self.id, self.name)

        return msg

    def printInfo(self):
        
        print(self)
        print("Elastic modulus: %.3f MPa" % self.E)
        print("Shear   modulus: %.3f MPa" % self.G)
        print("Density: %.3f g/mm^3" % self.D)
        print("Strain Limit: %.3f%%" % (self.max_strain * 100))

        print()

    def output(self):

        info = []
        info += [self.id]
        info += [self.name]
        info += [self.G]
        info += [self.E]
        info += [self.D]
        info += [self.max_strain]

        info = tuple(info)

        return info