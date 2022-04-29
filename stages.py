from parameters import MAXDOFS, EUCSPACEDIM
from messages import MSG_UNIMPLEMENTED
from utils import fixFloat
import numpy as np
import math
import sympy as sp
from screwVectors import swapOperator, Frame, trMatrix
from scipy.linalg import null_space
from DOFsystem import DOFSystem
from simulation import PropertyMatrix

class Stage(object):
    def __init__(self, id, material, center, volume):
        self.id = id
        self.center = fixFloat(np.asarray(center))
        self.volume = volume
        self.frame = None
        
        self.actuations = []

        self.twistPrms = np.zeros(MAXDOFS) # does not change, name or boundary conditions
        self.twistVars = np.zeros(MAXDOFS) # numbers that would change over operation

        self.material = material

    def getPrms(self):

        prms = []

        for val in self.twistPrms:
            if(isinstance(val, sp.core.symbol.Symbol)):
                prms += [val]
        
        return prms
    
    def setPrms(self, prms):
            
        for i in range(len(self.twistPrms)):
            val = self.twistPrms[i]
            if(isinstance(val, sp.core.symbol.Symbol)):
                tarVal = prms[val]
            else:
                tarVal = val
            
            self.twistVars[i] = tarVal
    
    def resetPrms(self):

        for i in range(len(self.twistPrms)):
            self.twistVars[i] = self.twistPrms[i]
    
    def addActuation(self, actuation):

        self.actuations += [actuation]

    def T0(self, tarFrame):

        TB = self.twistVars
        dB = tarFrame.origin - self.center
        DB = np.zeros(MAXDOFS)
        DB[EUCSPACEDIM:] = dB
        TrTB = trMatrix(TB)
        T0 = TB + np.matmul(trMatrix(TB), DB) - DB

        return T0

    @staticmethod
    def model(info, material):

        objType = info[0]
        for elemType in Stage.__subclasses__():
            if objType == elemType.name:
                modeled = elemType(info, material)

        return modeled

class RigidSolid(Stage):
    name = "RigidSolid"

    def __init__(self, info, material):
        super().__init__(info[1], material, info[3], info[4])

        self.actingWrench = np.zeros(MAXDOFS)

        self.isFixed = bool(info[5])

        self._getFrame()
        self._prmInit()

    def __repr__(self):

        msg = "Rigid stage(ID:%d, MAT:%s)" % (self.id, self.material.name)

        return msg

    def printInfo(self):
        
        print(self)
        print("Center: %s" % str(self.center))
        print("Volume: %s mm^3" % str(self.volume))
        print("Is fixed: %s" % str(self.isFixed))

        print()

    def _getFrame(self):

        # find system
        # [o, x, y, z]

        sys = np.identity(3)
        frame = Frame(self.center, sys)
        self.frame = frame

    def _prmInit(self):

        prms = ["t%d_s%s" % (i, self.id) for i in range(MAXDOFS)]
        self.twistPrms = np.asarray(sp.var(prms, real=True))
        self.twistVars = np.zeros(self.twistPrms.shape)
    
    def setWrench(self, wrench):

        self.actingWrench = wrench

    def setFixedEnd(self):
        
        if(self.isFixed):
            self.twistPrms = np.zeros(MAXDOFS).astype(object)

    def negWork(self):
        
        swapOp = swapOperator()
        work = np.matmul(np.matmul(self.twistVars, swapOp), self.actingWrench)
        
        return -1 * work