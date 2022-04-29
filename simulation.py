import numpy as np
from scipy.linalg import null_space
from parameters import EUCSPACEDIM, USE_TIMOSHENKO, MAXDOFS, EPSILON, EUCSPACEDIM
from subspace import Subspace
from screwVectors import decompose, screwVec

class PropertyMatrix(object):
    def __init__(self, material, flexure):
        
        self.material = material
        self.flexure = flexure

    def compliance(self):

        E = self.material.E # elastic modulus
        G = self.material.G # shear modulus
        A = self.flexure.A() # cross section area
        l = self.flexure.l() # length of element
        J = self.flexure.J() # torsion constant
        _, Iy, Iz = self.flexure.I() # area moments of inertia

        dx_fx = l / (E * A)
        dy_fy = (l ** 3) / (3 * E * Iz)
        dz_fz = (l ** 3) / (3 * E * Iy)

        ty_fz = dz_my = -1 * (l ** 2) / (2 * E * Iy)
        tz_fy = dy_mz = (l ** 2) / (2 * E * Iz)
        
        tx_mx = l / (G * J)
        ty_my = l / (E * Iy)
        tz_mz = l / (E * Iy)

        if(USE_TIMOSHENKO):
            dy_fy += l / (G * A)
            dz_fz += l / (G * A)
        
        C = [[    0,     0,     0, tx_mx,     0,     0],
             [    0,     0, ty_fz,     0, ty_my,     0],
             [    0, tz_fy,     0,     0,     0, tz_mz],
             [dx_fx,     0,     0,     0,     0,     0],
             [    0, dy_fy,     0,     0,     0, dy_mz],
             [    0,     0, dz_fz,     0, dz_my,     0]]

        compliance = np.asarray(C)
        
        return compliance

    def stiffness(self):
        
        compliance = self.compliance()
        stiffness = np.linalg.inv(compliance)

        return stiffness
