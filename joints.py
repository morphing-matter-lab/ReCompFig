import numpy as np
import math
from scipy.linalg import null_space
import sympy as sp

from DOFsystem import DOFSystem
from simulation import PropertyMatrix
from parameters import SUB_DIV_COUNT, MAXDOFS, EUCSPACEDIM, EPSILON
from messages import MSG_UNIMPLEMENTED
from screwVectors import screwTransform, screwVec, adjTransMatrix, findTrans, swapOperator, Frame, inverseAd, trMatrix
from utils import fixFloat, exhaustiveSubclasses
from stages import RigidSolid, Stage
from subspace import AxialSubspace, Subspace

class Joint(object):
    def __init__(self, id, adj):
        self.id = id
        self.adj = adj

class Flexure(object):
    name = "Flexure"
    isReal = True

    def __init__(self, id, adj, DOF, material, isSeg):

        self.id = id
        self.adj = adj
        self.freedom = None
        self.constraint = None
        self.DOF = None
        self.frameStart = None
        self.frameEnd = None

        self.base = None
        self.targ = None
        
        self.segments = []
        self.isSeg = isSeg

        self.twistPrms = np.zeros(MAXDOFS) # does not change, name or boundary conditions
        self.twistVars = np.zeros(MAXDOFS) # numbers that would change over operation
        
        self.material = material
        self.property = PropertyMatrix(material, self)
        
        if(isSeg):
            self._segInit()

    def _segInit(self):

        prms = ["t%d_j%s" % (i, self.id) for i in range(MAXDOFS)]
        self.twistPrms = np.asarray(sp.var(prms, real=True))
        self.twistVars = np.zeros(self.twistPrms.shape)

    def compliance(self):
        
        return self.property.compliance()

    def stiffness(self):
        
        return self.property.stiffness()

    def chain(self, startObj, endObj):

        self.base = startObj
        self.targ = endObj

    def getPrms(self):

        prms = []

        if(self.isSeg): # base case
            for val in self.twistPrms:
                if(isinstance(val, sp.core.symbol.Symbol)):
                    prms += [val]

        else: # recursion case
            for seg in self.segments:
                prms += seg.getPrms()
            
        return prms
    
    def setPrms(self, prms):

        if(self.isSeg): # base case
            for i in range(len(self.twistPrms)):
                val = self.twistPrms[i]
                if(isinstance(val, sp.core.symbol.Symbol)):
                    tarVal = prms[val]
                else: 
                    tarVal = val
                self.twistVars[i] = tarVal

        else: # recursion case
            for seg in self.segments:
                seg.setPrms(prms)
    
    def resetPrms(self):
        
        if(self.isSeg):
            for i in range(len(self.twistPrms)):
                self.twistVars[i] = self.twistPrms[i]
        else:
            for seg in self.segments:
                seg.resetPrms()
    
    def setVariables(self, valMap):
        
        for seg in self.segments:
            
            for i in range(len(seg.twistPrms)):
                if(isinstance(seg.twistPrms[i], sp.core.symbol.Symbol)):
                    seg.twistPrms[i] = valMap[seg.twistPrms[i]]
            
            seg.twistPrms = seg.twistPrms.astype(np.float)

    def collectTrans(self):
        
        cur  = self
        rooted = isinstance(cur, Stage)
        transSum = None
        while (not rooted):
            prev = cur.base

            # get Tn depending on base object type
            if(isinstance(prev, Stage)): T = prev.T0(cur.frameStart)
            else: T = prev.twistVars

            if(not isinstance(transSum, np.ndarray)):
                transSum = T + 0 # avoid aliasing by using add
            else:
                transSum += T
            
            cur = prev
            rooted = isinstance(cur, Stage)
        
        return transSum

    def strainEnergy(self):
        
        # recursion case: not segment -> iterate on segments
        if(not self.isSeg):
            result = 0
            for seg in self.segments:
                result += seg.strainEnergy()
            return result
        
        # base case: calculate strain energy within element
        # self twist parameters
        t = self.twistVars
        # swap operator
        swap = swapOperator()
        
        # find ad (self to world)
        transSum = self.collectTrans()
        transFrame = self.frameStart.transform(transSum)
        ad = transFrame.AdToWorld()
        adInv = inverseAd(ad)

        kElem = self.stiffness() # in elemental frame
        adElem2Frame = adjTransMatrix(Frame.spatial(), transFrame)
        adElem2FrameInv = inverseAd(adElem2Frame)
        k = np.matmul(np.matmul(adElem2Frame, kElem), adElem2FrameInv) # in current frame
        multiplied, components = t.T, [swap, ad, k, adInv, t]
        for i in range(len(components)):
            multiplied = np.matmul(multiplied, components[i])
        
        return .5 * multiplied

    def kinematicConstraints(self):
        
        T0 = self.base.T0(self.segments[0].frameStart)
        TsegSum = np.zeros(MAXDOFS)
        for seg in self.segments:
            TsegSum += seg.twistVars
        Tend = T0 + TsegSum

        Tr = trMatrix(Tend)

        dR = self.targ.center - self.frameEnd.origin
        DR = np.zeros(MAXDOFS)
        DR[EUCSPACEDIM:] = dR

        TR = self.targ.twistVars
        constraintEq = Tend + np.matmul(Tr, DR) - DR - TR
        
        return np.linalg.norm(constraintEq) # must equal to 0
    
    def isAllowed(self, consVecs):
        # assuming that consVecs are orthogonal
        
        concat = np.concatenate([self.constraint, consVecs])
        consSpace = Subspace(spans=consVecs)
        unionSpace = Subspace(spans=concat)
        
        return consSpace.spanRank == unionSpace.spanRank

    def dirPosAllowed(self, consVecs):

        dir, pos = True, True
        for vec in self.constraint:
            vec = vec / np.linalg.norm(vec)
            vec = vec.reshape((1, -1))
            dotProd = np.sum(vec * consVecs, axis=-1)
            projected = np.matmul(dotProd, consVecs)
            diff = vec - projected
            dirResidual = np.linalg.norm(diff[:EUCSPACEDIM])
            posResidual = np.linalg.norm(diff[EUCSPACEDIM:])
            if(dirResidual > EPSILON):
                dir = False
            if(posResidual > EPSILON):
                pos = False
        
        return dir, pos

    def isAssociated(self, axialSubspace):

        selfAS = self.axialSubspace()
        # check axial correlation
        axialCorr = False
        for axis in selfAS.axis:
            dotProd = np.sum(axis * axialSubspace.axis, axis=-1)
            isCorr = np.linalg.norm(dotProd) > EPSILON
            if(isCorr):
                axialCorr = True
                break
        
        if(not axialCorr): return False

        # check spatial correlation
        intersection = selfAS.intersect(axialSubspace)
        
        if(intersection == None): return False # not correlated at all
        if(axialSubspace.degree > axialSubspace.spanRank):
            spatialCorr = intersection != None
        else:
            spatialCorr = intersection.spanRank >= 1

        corr = axialCorr and spatialCorr
        
        return corr
    
    def screwSpace(self):

        return self.axialSubspace().outputScrewSpace()

    ############################################################################
    # inherited
    ############################################################################

    def A(self):
        # placeholder
        # returns the cross-section area of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def l(self):
        # placeholder
        # returns the length of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def J(self):
        # placeholder
        # returns the torsion constant of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def I(self):
        # placeholder
        # returns the area moments of inertia of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None, None, None

    def Ad(self):

        Ad = adjTransMatrix(self.frameEnd, self.frameStart)

        return Ad
    
    def Tr(self, refFrame):

        R = findTrans(self.frameEnd.system, refFrame.system)
        Tr =  np.zeros(MAXDOFS * 2, MAXDOFS * 2)
        Tr[0:3, 0:3] = R
        Tr[3:6, 3:6] = R
        Tr[6:9, 6:9] = R
        Tr[9:12, 9:12] = R
        
        return Tr

    def K(self):

        Ad = self.Ad()
        k = self.stiffness()
        K = np.zeros(MAXDOFS * 2, MAXDOFS * 2)
        K[:MAXDOFS, :MAXDOFS] = np.matmul(Ad.T, np.matmul(k, Ad))
        K[:MAXDOFS, MAXDOFS:] = -1 * np.matmul(Ad.T, k)
        K[MAXDOFS:, :MAXDOFS] = -1* np.matmul(k, Ad)
        K[MAXDOFS:, MAXDOFS:] = k

        return K

    def DB(self):
        # placeholder
        # returns the base displacement vector of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def DR(self):
        # placeholder
        # returns the target displacement vector of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def subDivide(self):
        # placeholder
        # returns a list of subdivided segments of the flexure element

        return [] # list of elements

    def axialSubspace(self):
        # placeholder
        # returns an AxialSubspace

        return 42

    ############################################################################
    # static
    ############################################################################

    @staticmethod
    def model(info, material):

        objType = info[0]
        for elemType in exhaustiveSubclasses(Flexure):
            if objType == elemType.name:
                modeled = elemType(info, material)

        return modeled

class Imaginary(Flexure):
    name = "Imaginary"
    isReal = False
    DOF = 6

    def __init__(self, info, material, isSeg=False):
        super().__init__(info[1], info[2], PlaceHolder.DOF, material, isSeg)
        self.freedom = np.identity(MAXDOFS)
        self.constraint = np.zeros((0, 0))

class WireCircular(Flexure):
    name = "WireCircular"
    DOF = 5

    def __init__(self, info, material, isSeg=False):
        super().__init__(info[1], info[2], WireCircular.DOF, material, isSeg)

        self.start = np.asarray(info[4]).astype(np.float)
        self.end = np.asarray(info[5]).astype(np.float)
        self.length = np.linalg.norm(self.end - self.start)
        self.radius = info[6]
        
        self._getSpaces()
        self._getFrames()

    def __repr__(self):

        msg = "Circular wire flexure(ID:%.2f, MAT:%s)" % (self.id, self.material.name)

        return msg

    def printInfo(self):
        
        print(self)
        print("Start: %s" % str(self.start))
        print("End: %s" % str(self.end))
        vec = self.end - self.start
        vec /= np.linalg.norm(vec)
        vec = fixFloat(vec)
        print("Vector: %s" % str(vec))
        print("Length: %.3f mm" % self.length)
        print("Radius: %.3f mm" % self.radius)
        print("Base: %s" % str(self.base))
        print("Targ: %s" % str(self.targ))

        print()

    def _getSpaces(self):
        
        refPt = fixFloat(self.start)
        vec = fixFloat(self.end - self.start)
        self.constraint = screwVec(vec, refPt).reshape((1, -1))
        self.freedom = DOFSystem.conToFree(self.constraint)

    def _getFrames(self):

        # find system
        # [o, x, y, z]

        xVec = (self.end - self.start) / self.length
        xVec = fixFloat(xVec)
        yzVecs = null_space(xVec.reshape(1, -1)).T
        yVec = yzVecs[0]
        zVec = np.cross(xVec, yVec)
        sys = np.stack([xVec, yVec, zVec])
        self.frameStart = Frame(self.start, sys)
        self.frameEnd = Frame(self.end, sys)

    ############################################################################
    # inherited
    ############################################################################

    def A(self):
        
        area = (self.radius ** 2) * math.pi

        return area

    def l(self):

        return self.length

    def J(self):
        # placeholder
        # returns the torsion constant of the flexure element

        J = (self.radius ** 4) * math.pi / 2

        return J

    def I(self):
        # placeholder
        # returns the area moments of inertia of the flexure element

        Ix = None
        Iy = Iz = (self.radius ** 4) * math.pi / 4

        return Ix, Iy, Iz

    def DB(self):
        # placeholder
        # returns the base displacement vector of the flexure element

        vec = self.frameStart[0] - self.base.center

        return vec

    def DR(self):
        # placeholder
        # returns the target displacement vector of the flexure element

        vec = self.targ.center - self.frameEnd[0]

        return vec

    def subDivide(self, segCount=SUB_DIV_COUNT):

        # generate list of nodes
        ids = [self.adj[0]]
        for i in range(segCount - 1):
            id = self.id + ((i + 1) * 0.01)
            ids += [id]
        ids += [self.adj[1]]

        # generate new elements
        vec = self.end - self.start
        self.segments = []
        for i in range(segCount):
            startPt = self.start + (vec * (i / segCount))
            endPt = self.start + (vec * ((i + 1) / segCount))

            info = []
            info += [WireCircular.name]
            info += [self.id + ((i + 1) * 0.01)] # id
            info += [(ids[i], ids[i + 1])] # adj, does not matter
            info += [self.material.id] # material
            info += [startPt] # start point
            info += [endPt] # end point
            info += [self.radius]

            newElem = WireCircular(info, self.material, isSeg=True)
            self.segments += [newElem]
        
        # chain elements
        for i in range(len(self.segments)):
            if (i == 0): base = self.base
            else: base = self.segments[i - 1]
            
            if (i == (len(self.segments) - 1)): targ = self.targ
            else: targ = self.segments[i + 1]

            self.segments[i].chain(base, targ)

        return self.segments

    def axialSubspace(self):

        axis = (self.end - self.start) / self.length
        refPt = self.start
        spanVecs = axis.reshape((-1, EUCSPACEDIM))
        
        return AxialSubspace(axis, refPt, spanVecs)

class Cable(WireCircular):
    name = "Cable"
    DOF = 5

    def __init__(self, info, material, isSeg=False):
        super().__init__(info, material, isSeg=isSeg)

    def __repr__(self):

        msg = "Cable flexure(ID:%.2f, MAT:%s)" % (self.id, self.material.name)

        return msg

class Sensor(WireCircular):
    name = "Sensor"
    DOF = 6

    def __init__(self, info, material, isSeg=False):
        super().__init__(info, material, isSeg=isSeg)

    def __repr__(self):

        msg = "Sensor flexure(ID:%.2f, MAT:%s)" % (self.id, self.material.name)

        return msg

class PlaceHolder(Flexure):
    name = "PlaceHolder"
    isReal = False
    DOF = 6

    def __init__(self, info, material, isSeg=False):
        super().__init__(info[1], info[2], PlaceHolder.DOF, material, isSeg)
        self.center = np.asarray(info[4]).astype(np.float)

        self.freedom = np.identity(MAXDOFS)
        self.constraint = np.zeros((0, 0))

    def A(self):
        # placeholder
        # returns the cross-section area of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def l(self):
        # placeholder
        # returns the length of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def J(self):
        # placeholder
        # returns the torsion constant of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def I(self):
        # placeholder
        # returns the area moments of inertia of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None, None, None

    def Ad(self):

        Ad = adjTransMatrix(self.frameEnd, self.frameStart)

        return Ad
    
    def Tr(self, refFrame):

        R = findTrans(self.frameEnd.system, refFrame.system)
        Tr =  np.zeros(MAXDOFS * 2, MAXDOFS * 2)
        Tr[0:3, 0:3] = R
        Tr[3:6, 3:6] = R
        Tr[6:9, 6:9] = R
        Tr[9:12, 9:12] = R
        
        return Tr

    def K(self):

        Ad = self.Ad()
        k = self.stiffness()
        K = np.zeros(MAXDOFS * 2, MAXDOFS * 2)
        K[:MAXDOFS, :MAXDOFS] = np.matmul(Ad.T, np.matmul(k, Ad))
        K[:MAXDOFS, MAXDOFS:] = -1 * np.matmul(Ad.T, k)
        K[MAXDOFS:, :MAXDOFS] = -1* np.matmul(k, Ad)
        K[MAXDOFS:, MAXDOFS:] = k

        return K

    def DB(self):
        # placeholder
        # returns the base displacement vector of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def DR(self):
        # placeholder
        # returns the target displacement vector of the flexure element

        assert False, MSG_UNIMPLEMENTED

        return None

    def subDivide(self):
        # placeholder
        # returns a list of subdivided segments of the flexure element

        return None # list of elements
