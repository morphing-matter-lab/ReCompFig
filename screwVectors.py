import numpy as np
from scipy.optimize.optimize import vecnorm
from parameters import MAXDOFS, EUCSPACEDIM, EPSILON
from utils import fixFloat, p2t
from sympy import symbols, sin, cos
import sympy as sp

class Frame(object):
    def __init__(self, origin, sys):
        self.origin = origin
        self.system = sys
        self.x = sys[0]
        self.y = sys[1]
        self.z = sys[2]

    def __repr__(self):

        msg = "O: %s X: %s Y: %s Z: %s" % (str(self.origin), str(self.x), str(self.y), str(self.z))

        return msg

    def transform(self, screw):

        # transform vectors
        origin = screwTransform(self.origin, screw)
        sys = np.stack([screwTransform(vec, screw, isPt=False) for vec in self.system])
        newFrame = Frame(origin, sys)

        return newFrame

    def output(self):
        
        o = p2t(self.origin)
        x = p2t(self.x)
        y = p2t(self.y)
        z = p2t(self.z)

        return tuple([o, x, y, z])

    def AdToWorld(self):

        Ad = adjTransMatrix(self, Frame.spatial())

        return Ad

    @staticmethod
    def spatial():
        return Frame(np.zeros(EUCSPACEDIM), np.identity(EUCSPACEDIM))

    @staticmethod
    def fromList(inp):
        o = np.asarray(inp[0])
        sys = np.asarray(inp[1:])

        return Frame(o, sys)
        
def swapOperator():

    swap = np.eye(MAXDOFS)
    swap = np.concatenate((swap[EUCSPACEDIM:, :], 
                           swap[:EUCSPACEDIM, :]), axis=0)
    
    return swap

def screwVec(axis=[0, 0, 0], ref=[0, 0, 0], trans=[0, 0, 0]):

    if(not isinstance(axis, np.ndarray)):
        axis = np.array(axis)
    if(not isinstance(ref, np.ndarray)):
        ref = np.array(ref)
    if(not isinstance(trans, np.ndarray)):
        trans = np.array(trans)

    s = np.cross(ref, axis) + trans
    
    result = np.concatenate((axis, s))

    return result

def lineToScrewVec(motion):

    # chr,      2 points
    motionType, p0, p1 = motion
    p0, p1 = np.asarray(p0), np.asarray(p1)
    vec = p1 - p0

    if(motionType == 'r' or motionType == 'w'): # rotation & wrench
        sv = screwVec(axis=vec, ref=p0)
    elif(motionType == 't'): # translation
        sv = screwVec(trans=vec)
    else:
        assert False, "function not implemented"

    sv = fixFloat(sv)
    return sv

def skewSymmetricMatrix(vec):

    x, y, z = vec
    m = np.asarray([[ 0, -z,  y],
                    [ z,  0, -x],
                    [-y,  x,  0]])

    return m

def getAxis(vec, unitize=False):

    vecNorm = np.linalg.norm(vec[:EUCSPACEDIM])

    if(vecNorm > EPSILON):
        if(unitize):
            return vec[:EUCSPACEDIM] / vecNorm
        else:
            return vec[:EUCSPACEDIM]
    else:
        if(unitize):
            return vec[EUCSPACEDIM:] / np.linalg.norm(vec[EUCSPACEDIM:])
        else:
            return vec[EUCSPACEDIM:]

def decompose(vec):

    # axis
    rotMag =   ((vec[0] ** 2) + (vec[1] ** 2) + (vec[2] ** 2)) ** .5
    transMag = ((vec[3] ** 2) + (vec[4] ** 2) + (vec[5] ** 2)) ** .5
    
    if(rotMag != 0): # rotation
        axis, mag = vec[:3], rotMag
        axisUnit = axis / mag
        transLen = np.dot(axis, vec[3:])
        crossProd = vec[3:] - transLen * axis
        refPt = np.cross(axisUnit, crossProd / mag)
    elif(transMag != 0): # translation
        axis, mag = vec[3:], 0
        axisUnit = axis / transMag
        transLen = transMag
        refPt = np.zeros(EUCSPACEDIM)
    else: # all 0
        axis, mag = np.zeros((EUCSPACEDIM, EUCSPACEDIM)), 0
        axisUnit = np.zeros(EUCSPACEDIM)
        transLen = 0
        refPt = np.zeros(EUCSPACEDIM)
        
    return axisUnit, mag, transLen, refPt

def screwTransform(vec, screwVec, isPt=True):
    
    if(not isinstance(vec, np.ndarray)):
        vec = np.asarray(vec)

    if(vec.shape[-1] == EUCSPACEDIM): # point vector mode
        axis, rotMag, transLen, refPt = decompose(screwVec)
        
        R = screwRot(axis, rotMag)
        if(isPt): # point
            transformed = refPt + np.matmul(R, vec - refPt) + transLen * axis
        else: # vec
            transformed = np.matmul(R, vec)

    else: # screw vector mode
        vecAxis, vecMagnitude, vecTransLen, vecRefPt = decompose(vec)
        axisNew = screwTransform(vecAxis, screwVec, isPt=False) # unit
        refPtNew = screwTransform(vecRefPt, screwVec, isPt=True)
        transformed = screwVec(axisNew * vecMagnitude, 
                               refPtNew, 
                               vecTransLen * axisNew)
    
    return transformed

def screwTransMat(d, R):
    D = skewSymmetricMatrix(d)
    m = np.zeros((MAXDOFS, MAXDOFS)).astype(d.dtype)
    m[:3, :3] = R
    m[3:, :3] = np.matmul(D, R)
    m[3:, 3:] = R

    return m

def findTrans(source, target):

    # Ax = B
    # A = Bx^-1
    # A: to find, B: target X: source

    trans = np.matmul(target, source.T)

    return trans

def adjTransMatrix(oldSys, newSys):

    oldPt = oldSys.origin
    oldFrame = oldSys.system

    newPt = newSys.origin
    newFrame = newSys.system

    R = findTrans(oldFrame, newFrame)
    d = oldPt - newPt
    d = np.matmul(newFrame, d)
    Ad = screwTransMat(d, R) # negative no rot positive w/rot

    return Ad

def screwRot(inputVec, rotMag=None):
    
    if(len(inputVec) == MAXDOFS): # screw mode, decompose
        axis, rotMag, _, _ = decompose(inputVec)
    elif(len(inputVec) == EUCSPACEDIM and rotMag != None): # axis and rotMag
        axis = inputVec
    K = skewSymmetricMatrix(axis)
    # watch out for types casting here
    R = np.identity(3) +\
        np.sin(rotMag) * K +\
        (1 - np.cos(rotMag)) * np.matmul(K, K)
        
    return R

def inverseAd(Ad):
    R = Ad[:EUCSPACEDIM, :EUCSPACEDIM]
    DR = Ad[EUCSPACEDIM:, :EUCSPACEDIM]
    D = np.matmul(DR, R.T) # r is rotation matrix -> transpose = inverse

    adInv = np.zeros((MAXDOFS, MAXDOFS))
    adInv[:EUCSPACEDIM, :EUCSPACEDIM] = R.T
    adInv[EUCSPACEDIM:, :EUCSPACEDIM] = -1 * np.matmul(R.T, D)
    adInv[EUCSPACEDIM:, EUCSPACEDIM:] = R.T

    return adInv

def trMatrix(screw):

    R = screwRot(screw)
    Tr = np.zeros((MAXDOFS, MAXDOFS))
    Tr[EUCSPACEDIM:, EUCSPACEDIM:] = R
    
    return Tr

def screwVecLen(vec):

    axis = vec[:EUCSPACEDIM]
    axisLen = np.linalg.norm(axis)
    if(axisLen > EPSILON): # rotation axis exists
        vecLen  = axisLen
    else:
        vecLen = np.linalg.norm(vec[EUCSPACEDIM:])
    
    return vecLen

def normScrew(vec):

    vecLen = screwVecLen(vec)
    
    if(vecLen <= EPSILON):
        return np.zeros(MAXDOFS)
    else:
        return vec / vecLen
