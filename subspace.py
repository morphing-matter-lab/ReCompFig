import numpy as np
from numpy.random import normal
import sympy as sp
import copy
from scipy.linalg import null_space
from itertools import combinations

from sympy.core.function import Subs

from screwVectors import screwVec, normScrew
from parameters import EPSILON, EUCSPACEDIM, SANITYCHECK, EPSILON2
from utils import fixFloat, solve

class Subspace(object):
    def __init__(self, normals=None, refPt=None, parameters=None, spans=None, dim=None):
        # null space as row vector (N, M), N: number of null spaces, M: dimension
        # ref pt is a point on the subspace (D): dimensions
        # parameters are the variables(symbols) of the subspace (M), M: number of parameters

        self.dim = dim

        self.normalSpace = None
        self.spanSpace = None
        self.refPt = None

        self.parameters = None
        self.isPoint = None
        self.isFull = None
        self.normalRank = 0
        self.spanRank = 0

        self.subsetSpaces = [] # for multimodal systems

        self._initialize(normals, spans, refPt, parameters)
        
    def __eq__(self, other):
        
        # atrribute check
        if(not isinstance(other, Subspace)): return False
        if(self.dim != other.dim): return False
        if(self.normalRank != other.normalRank): return False
        if(self.spanRank != other.spanRank): return False

        # check points (shortcut)
        if(self.isPoint):
            dist = abs(np.linalg.norm(self.refPt - other.refPt))
            return dist <= EPSILON
        
        # check if multually share position
        selfPtOnOther = other.isPtOnSubspace(self.refPt)
        otherPtOnSelf = self.isPtOnSubspace(other.refPt)
        if((not selfPtOnOther) or (not otherPtOnSelf)): return False

        # check if normal/span sapces are mutually complementary
        complementary = Subspace.isComplementary(self.normalSpace, other.spanSpace) and\
                        Subspace.isComplementary(self.spanSpace, other.normalSpace)
        if(not complementary): return False

        return True
        
    def __hash__(self):

        infoList = [self.dim, self.normalRank, self.spanRank]
        infoList += list(self.normalSpace.reshape(-1)) + list(self.spanSpace.reshape(-1))
        return hash(infoList)
    
    def __repr__(self):

        msg = "SubSpace%dD(Normals:%d, Spans:%d)" %\
              (self.dim, self.normalRank, self.spanRank)

        return msg

    def _initialize(self, normals, spans, refPt, parameters):
        
        normals, spans, refPt, parameters = self._formatInput(normals, spans, refPt, parameters)
        
        if(self.dim != None):
            self._initFullSpace()
        if(isinstance(normals, np.ndarray)): # construct using normals
            self.dim = normals.shape[-1]
            self._initByNormal(normals)
        elif(isinstance(spans, np.ndarray)): # construct using spans
            self.dim = spans.shape[-1]
            self._initBySpan(spans)
        
        self._initPoint(refPt, parameters)
        
        self._updateInfo()

    def _formatInput(self, normals, spans, refPt, parameters):

        # normals
        if(isinstance(normals, np.ndarray) and normals.size != 0):
            normals = normals.reshape((-1, normals.shape[-1]))
        elif(isinstance(normals, list) or isinstance(normals, tuple)):
            normals = np.asarray(normals)
        else:
            normals = None
        
        # spans 
        if(isinstance(spans, np.ndarray) and spans.size != 0):
            spans = spans.reshape((-1, spans.shape[-1]))
        elif(isinstance(spans, list) or isinstance(spans, tuple)):
            spans = np.asarray(spans)
        else:
            spans = None
        
        # reference point
        if(isinstance(refPt, np.ndarray)):
            refPt = refPt.reshape(-1)
        elif(isinstance(refPt, list) or isinstance(refPt, tuple)):
            refPt = np.asarray(refPt)
        else:
            refPt = None

        # parameters
        if(isinstance(parameters, np.ndarray)):
            parameters = parameters.reshape(-1)
        elif(isinstance(parameters, list) or isinstance(parameters, tuple)):
            parameters = np.asarray(parameters)
        else:
            parameters = None

        return normals, spans, refPt, parameters

    def _initFullSpace(self):

        self.normalSpace = np.asarray([[]])
        self.spanSpace = np.eye(self.dim)
        
    def _initByNormal(self, normals):
        
        space = Subspace(dim=self.dim)
        for vec in normals: 
            space.expandNormal(vec)

        self._getDataFrom(space)

    def _initBySpan(self, spans):

        space = Subspace(refPt=np.zeros(self.dim))
        for vec in spans: space.expandSpan(vec)
        
        self._getDataFrom(space)

    def _initPoint(self, refPt, parameters):

        if(self.dim == None):
            # dimension unknown, figure out with refPt and parameters
            if(isinstance(refPt, np.ndarray)): 
                self.dim = refPt.shape[-1]
            elif(isinstance(parameters, np.ndarray)): 
                self.dim = parameters.shape[-1]

        # reference point
        if(not isinstance(refPt, np.ndarray)):
            self.refPt = np.zeros(self.dim)
        else:
             self.refPt = refPt
             
        # parameters
        if(not isinstance(parameters, np.ndarray)):
            names = ' '.join([chr(ord("a") + i) for i in range(self.dim)])
            self.parameters = np.asarray(sp.symbols(names))
        else:
            self.parameters = parameters
        
        # span
        if(not isinstance(self.spanSpace, np.ndarray)):
            self.spanSpace = np.asarray([[]])
        
        # normals
        if(not isinstance(self.normalSpace, np.ndarray)):
            self.normalSpace = np.eye(self.dim)
        
    def _getDataFrom(self, other):

        self.dim = other.dim

        self.normalSpace = other.normalSpace
        self.spanSpace = other.spanSpace
        self.parameters = other.parameters
        self.refPt = other.refPt

        self._updateInfo()

    def _sanityCheck(self):

        # check normals and spans
        assert Subspace.isComplementary(self.normalSpace, self.spanSpace),\
               "normal and span space are not orthogonal"
        
        # check ref point on space
        assert self.isPtOnSubspace(self.refPt), "refPt not on plane"

        # check if subspace has parameters
        assert self.parameters.size == self.dim, "parameters not right"

        # check all vectors are unitized
        bases = self.getBases()
        length = np.linalg.norm(bases, axis=-1)
        deviation = abs(length) - 1
        assert np.all(deviation <= EPSILON), "bases are not unitized"

        # check all vectors are orthogonal
        for i in range(self.dim):
            for j in range(self.dim):
                if(i < j):
                    dotProd = np.dot(bases[i], bases[j])
                    assert abs(dotProd) <= EPSILON, "bases not orthogonal"

    def _updateInfo(self):

        self.spanRank = int(self.spanSpace.size / self.dim)
        self.normalRank = int(self.normalSpace.size / self.dim)
        self.isPoint = self.spanRank == 0
        self.isFull = self.normalRank == 0

        if(SANITYCHECK): self._sanityCheck()

    def _intersect(self, other):
        # note: other allways has a lower or equal solution space rank than self
        
        if(self.isParallel(other)): # two subspaces are parallel
            return None
        elif(self.includesSubspace(other)):
                return other
        
        # find intersection subspace by expanding normals
        normals = self.normalSpace
        newNormals = []
        for vec in other.normalSpace:
            # expand normal
            expanded = Subspace.expand(normals, vec)
            if(expanded.shape == normals.shape): # nothing happened
                continue
            newNormals += [expanded[-1]]
            normals = expanded
        
        newNormals = np.stack(newNormals)
        
        if(self.isPtOnSubspace(self.refPt) and other.isPtOnSubspace(self.refPt)):
            refPt = self.refPt
        else:
            symbs = sp.symbols(' '.join(["t%d" % i for i in range(len(newNormals))]))
            symbs = np.asarray(symbs).reshape(-1)
            correction = np.sum(newNormals * symbs.reshape((-1, 1)), axis=0)
            corrected = self.refPt + correction
            pointerVec = other.refPt - corrected
            otherNormalDot = np.dot(other.normalSpace, pointerVec)
            eqs = []
            for dotProd in otherNormalDot:
                if(len(dotProd.free_symbols) == 0):
                    if(abs(dotProd) <= EPSILON): continue
                    else: return None
                else:
                    eqs += [sp.Eq(dotProd, 0)]

            factor = solve(eqs, symbs)
            refPt = self.refPt + np.sum(newNormals * factor.reshape(-1, 1), axis=0)
            
            if(SANITYCHECK):
                assert (self.isPtOnSubspace(refPt) and\
                        other.isPtOnSubspace(refPt)), "implementation error"
            
        # constrcut new subspace
        intersection = Subspace(normals=normals, 
                                refPt=refPt, 
                                parameters=self.parameters)
        
        return intersection
    
    def _isParallelVec(self, span):

        if(self.isFull):
            return False
        
        for vec in span:
            dotProd = np.dot(self.normalSpace, vec)
            if(np.any(abs(dotProd) > EPSILON)):
                return False

        return True

    def _isParallelSubspace(self, other):

        pointer = other.refPt - self.refPt
        # get all span vectors
        catList = []
        if(self.spanRank > 0): catList += [self.spanSpace]
        if(other.spanRank > 0): catList += [other.spanSpace]
        if(len(catList) == 0): # shortcut: both are points, evaluate distance
            return np.linalg.norm(pointer) > EPSILON
        spans = np.concatenate(catList, axis=0)

        spansOrtho = np.asarray([[]])
        for vec in spans: spansOrtho = Subspace.expand(spansOrtho, vec)

        correction = Subspace.project(pointer, spansOrtho)
        residual = pointer - correction
        isParallel = np.any(abs(residual) > EPSILON) # there are components that cannot be eliminated

        return isParallel

    def _cpPoint(self, pt):

        if(self.isFull): return pt # full space special case
        
        vec = self.refPt - pt
        scale = np.dot(self.normalSpace, vec).reshape((-1, 1))
        corrections = self.normalSpace * scale
        correction = np.sum(corrections, axis=0)
        
        closestPt = pt + correction

        if(SANITYCHECK):
            assert self.isPtOnSubspace(closestPt), "implementation error"

        return closestPt

    def _cpSubspace(self, other):

        if(self.isFull): # handle edge cases
            if(other.isPoint): return other.refPt
            else: return other
        
        if(not self.isParallel(other)): # two subspaces would intersect
            return self.intersect(other)
        else: # tow subspaces are parallel
            # step 1: pull other to self
            cloned = Subspace.clone(other)
            pointer = self.refPt - other.refPt
            # eliminate the "normal" part from both subspaces
            spansCollected = self.spanSpace
            for vec in other.spanSpace:
                spansCollected = Subspace.expand(spansCollected, vec)
            
            if(spansCollected.size != 0):
                spansProj = Subspace.project(pointer, spansCollected)
            else:
                spansProj = np.zeros(pointer.shape)

            correction = pointer - spansProj
            cloned.move(correction)

            # find intersection
            closest = self.intersect(cloned)

            return closest

    def move(self, vec):

        self.refPt = self.refPt + vec

    def expandSpan(self, vec):
        
        # check if vectorhas zero length
        if(np.linalg.norm(vec) <= EPSILON):
            return # no action needed

        expanded = Subspace.expand(self.spanSpace, vec)
        if(self.spanSpace.size == expanded.size): # no changes, skip
            pass
        
        isFull = expanded.size / self.dim == self.dim
        if(isFull): normals = np.asarray([[]])
        else: normals = null_space(expanded).T
        
        self.spanSpace = expanded
        self.normalSpace = normals
        self._updateInfo()

    def expandNormal(self, vec):

        expanded = Subspace.expand(self.normalSpace, vec)
        if(self.normalSpace.shape == expanded.shape): # no changes, skip
            pass
        
        isFull = expanded.size / self.dim == self.dim
        if(isFull): spans = np.asarray([[]])
        else: spans = null_space(expanded).T

        self.spanSpace = spans
        self.normalSpace = expanded
        self._updateInfo()

    def removeSpan(self, vec):

        removed = Subspace.remove(self.spanSpace, vec)
        if(self.spanSpace.shape == removed.shape):
            pass

        isZero = removed.size == 0
        if(isZero): normals = np.eye(self.dim)
        else: normals = null_space(removed).T

        self.spanSpace = removed
        self.normalSpace = normals
        self._updateInfo()

    def removeNormal(self, vec):
        
        removed = Subspace.remove(self.normalSpace, vec)
        if(self.normalSpace.shape == removed.shape):
            pass

        isZero = removed.size == 0
        if(isZero): spans = np.eye(self.dim)
        else: spans = null_space(removed).T

        self.spanSpace = spans
        self.normalSpace = removed
        self._updateInfo()

    def printInfo(self, fullReport=False, printSpan=False, printPrm=False, \
                  printNormal=False, printRefPt=False, printTitle=False):

        print("Subspace entity (Dim: %d, Normals: %d, Spans: %d)" %\
              (self.dim, self.normalRank, self.spanRank))
        
        if(fullReport or printPrm):
            print("Parameters:")
            print(self.parameters)
        if(fullReport or printRefPt):
            print("point on subspace:")
            print(self.refPt)
        if(fullReport or printSpan):
            print("Spans:")
            print(self.spanSpace)
        if(fullReport or printNormal):
            print("Normals:")
            print(self.normalSpace)
        
        print()

    def includesSubspace(self, other):

        if(self.isFull):
            return True
        elif(self.isPoint and other.isPoint):
            return np.linalg.norm(self.refPt - other.refPt) <= EPSILON
        # check if span is included
        for vec in other.spanSpace:
            dotProd = np.dot(self.normalSpace, vec)
            if(np.any(abs(dotProd) > EPSILON)): return False
        
        # check if other conincide with self
        isOn = self.eval(other.refPt) <= EPSILON

        return isOn

    def isSubspaceOf(self, other):

        return other.includesSubspace(self)

    def isParallel(self, other):
        # parallel: two objects do not intersect

        if(isinstance(other, Subspace)):
            span = other.spanSpace
            vecInput = False
        else: # a vector
            span = other.reshape((-1, self.dim))
            vecInput = True
        
        if(vecInput):
            return self._isParallelVec(span)
        else:
            return self._isParallelSubspace(other)

    def isPtOnSubspace(self, pt):
        
        if(self.normalRank != 0):
            dist = abs(self.eval(pt))
        else:
            dist = np.linalg.norm(self.refPt - pt)
        
        onPlane = dist <= EPSILON

        return onPlane

    def closestPoint(self, other):
        
        if(isinstance(other, np.ndarray)):
            return self._cpPoint(other)
        elif(isinstance(other, Subspace)):
            return self._cpSubspace(other)

    def minDistTo(self, other):
        
        if(isinstance(other, np.ndarray)):
            return abs(self.eval(other))
        elif(isinstance(other, Subspace)):
            cp1 = self.closestPoint(other)
            cp2 = other.closestPoint(self)
            if(cp1 == cp2): # intersection found
                return 0
            else:
                return abs(cp1.eval(cp2.refPt))
        
    def eval(self, params):
        
        dotProd = np.dot(self.normalSpace, params)
        evalDist = np.linalg.norm((dotProd.reshape(-1)))
        offSetDotProd = np.dot(self.normalSpace, self.refPt)
        offset = np.linalg.norm((offSetDotProd.reshape(-1)))

        result = evalDist - offset

        return result
    
    def eqs(self):

        equations = []
        for normal in self.normalSpace:
            lhs = np.dot(normal, self.parameters)
            rhs = np.dot(normal, self.refPt)
            eq = sp.Eq(lhs, rhs)
            equations += [eq]
        
        return equations

    def intersect(self, other):
        
        # make sure self has higher span rank than other
        if(self.spanRank < other.spanRank):
            return other.intersect(self)

        # cases dispathcer
        if(self.isFull): # self is full rank, simply return the other
            return other
        elif(self == other): # special case, self and other are identical
            
            return self
        elif(self.isPoint and other.isPoint): # both are points
            dist = abs(np.linalg.norm(self.refPt - other.refPt))
            if(dist <= EPSILON): return self
        elif((not self.isPoint) and other.isPoint):
            if(self.isPtOnSubspace(other.refPt)): return other
        else:
            return self._intersect(other)
        
        # float point edge case
        if(self.minDistTo(other) <= EPSILON2):
            # use cp
            return other.closestPoint(self)

        return None
    
    def output(self, tVecs, mode='r'):

        result = {}
        twist = np.matmul(tVecs.T, self.refPt)
        twist = normScrew(twist)
        axis, normal = twist[:EUCSPACEDIM], twist[EUCSPACEDIM:]
        ptOnAxis = fixFloat(np.cross(axis, normal))
        result["refPt"] = tuple(ptOnAxis)
        
        spanVecs = []
        if(self.spanRank > 0):
            # get all possible combincations of span vectors
            combs = []
            for i in range(self.spanRank):
                combs += list(combinations(self.spanSpace, i + 1))
            
            # gather all compound vectors to avoid "singularities" when 
            # projecting high dimensional spaces into 3D
            allSpanVecs = []
            for comb in combs:
                summed = np.sum(np.asarray(comb), axis=0)
                allSpanVecs += [summed]
            
            # project vectors into 3D
            for spanVec in allSpanVecs:
                twist = np.matmul(tVecs.T, spanVec)
                normal = twist[EUCSPACEDIM:]
                
                if(mode == 'r'):
                    vecIn3D = np.cross(axis, normal)
                elif(mode == 't'):
                    vecIn3D = normal
                else:
                    assert False, "function not implemented"
                # normalize
                vecNorm = np.linalg.norm(vecIn3D)
                if(vecNorm == 0): continue
                vecUnit = vecIn3D / vecNorm
                vec = fixFloat(vecUnit)
                vec = vec / np.linalg.norm(vec)
                spanVecs += [vec]
            result["spanVecs"] = tuple([tuple(spanVec) for spanVec in spanVecs])
        else:
            result["spanVecs"] = tuple()
        
        return result

    def outputRefPtSpan(self, spanOnly=False):

        result = [] if spanOnly else [self.refPt]

        for span in self.spanSpace:
            if (len(span) != self.dim):
                continue
            result += [span]
        
        result = np.stack(result)

        return result

    def getBases(self):

        if(self.isFull):
            bases = self.spanSpace
        elif(self.isPoint):
            bases = self.normalSpace
        else:
            bases = np.concatenate((self.spanSpace, self.normalSpace), axis=0)

        return bases

    def copy(self):

        if(self.spanSpace.size != 0):
            new = Subspace(spans=np.copy(self.spanSpace), refPt=np.copy(self.refPt))
        else:
            new = Subspace(normals=np.copy(self.normalSpace), refPt=np.copy(self.refPt))
        
        return new

    @staticmethod
    def project(vec, space):

        spaceNorm = np.linalg.norm(space, axis=-1).reshape((-1, 1))
        spaceUnit = space / spaceNorm
        dotProd = np.dot(spaceUnit, vec).reshape((-1, 1))
        vecs = spaceUnit * dotProd
        proj = np.sum(vecs, axis=0)

        return proj

    @staticmethod
    def expand(space, vec):

        vecNorm = np.linalg.norm(vec)
        if(abs(np.linalg.norm(vec)) <= EPSILON): 
            return space
        if(space.size == 0):
            return vec.reshape((1, -1))
        
        vecNormed = vec / vecNorm
        dotProd = np.dot(space, vecNormed)
        
        if(np.all(abs(dotProd) <= EPSILON)): # perpendicular to all
            catList = (space, vec.reshape((1, -1)))
            return np.concatenate(catList, axis=0)
        else:
            correction = Subspace.project(vecNormed, space)
            residual = vecNormed - correction
            if(np.all(abs(residual) <= EPSILON)): return space

            newVec = residual / np.linalg.norm(residual)
            catList = (space, newVec.reshape((1, -1)))
            
            return np.concatenate(catList, axis=0)
            
    @staticmethod
    def remove(space, vec):

        if(abs(np.linalg.norm(vec)) <= EPSILON): 
            return space
        if(space.size == 0):
            return space

        dotProd = np.dot(space, vec)
        if(np.all(abs(dotProd) <= EPSILON)): # perpendicular to all
            return space
        else:
            normalSpace = null_space(space).T
            normalExpanded = Subspace.expand(normalSpace, vec)
            if(normalExpanded.size == normalSpace.size): return space
            removedSpace = null_space(normalExpanded).T

            return removedSpace

    @staticmethod
    def isComplementary(space1, space2):
        
        if(space1.size == 0):
            rank = np.linalg.matrix_rank(space2)
            return rank == space2.shape[-1]
        elif(space2.size == 0):
            rank = np.linalg.matrix_rank(space1)
            return rank == space1.shape[-1]
        else:
            for vec in space2:
                dotProd = np.dot(space1, vec)
                if(np.any(abs(dotProd) > EPSILON)): return False
            
            s1NullRank = null_space(space1).size / space1.shape[-1]
            s2NullRank = null_space(space2).size / space2.shape[-1]
            rank = s1NullRank + s2NullRank
            isFullRank = rank == space1.shape[-1]

            return isFullRank

    @staticmethod
    def cullDuplicates(subspaces):

        cleaned = []
        for item in subspaces:
            if(not isinstance(item, Subspace)):
                cleaned += [item]
            elif(item not in cleaned):
                cleaned += [item]
        
        return subspaces

    @staticmethod
    def solve(subspaces, printInfo=SANITYCHECK):

        dim = subspaces[0].dim
        prms = subspaces[0].parameters

        # check if planes share an intersection
        # method: solve cummulative intersection
        solutionSpace = Subspace(dim=dim, parameters=prms)
        if(printInfo): print(solutionSpace)
        for subspace in subspaces:
            if(printInfo): print("\t", subspace)
            intersection = solutionSpace.intersect(subspace)
            if(printInfo): print(intersection)
            if(intersection == None):
                return None
            solutionSpace = intersection
        if(printInfo): print()
        
        return solutionSpace

    @staticmethod
    def clone(source):

        return copy.deepcopy(source)

    @staticmethod
    def checkParallel(first, second):
        if(isinstance(first, Subspace)):
            return first.isParallel(second)
        elif(isinstance(second, Subspace)):
            return second.isParallel(first)
        else:
            return Subspace(spans=first).isParallel(second)

    @staticmethod
    def difference(first, second, returnSubspace=True):

        if(not isinstance(first, Subspace) and isinstance(first, np.ndarray)):
            first = Subspace(spans=first)
        if(not isinstance(second, Subspace) and isinstance(second, np.ndarray)):
            first = Subspace(spans=second)
        
        # not(union(not A, B))
        diffInverse = np.concatenate([first.normalSpace, second.spanSpace])
        diffSpace = Subspace(normals=diffInverse)
        
        if(returnSubspace):
            return diffSpace
        else:
            return diffSpace.spanSpace

class AxialSubspace(Subspace):
    def __init__(self, axis, refPt, spanVecs, isPrinciple=True):

        self.axis = None
        self.degree = None

        self.isTrans = None
        self.isPrinciple = isPrinciple

        self._initAxialSubspace(axis, refPt, spanVecs)

    def __eq__(self, other):

        if(not isinstance(other, AxialSubspace)): return False

        subspaceSame = super().__eq__(other)
        if(not subspaceSame): return False

        # check axis
        if(self.isFull):
            return other.isFull
        else:
            axisNull = null_space(self.axis).T
            return Subspace.isComplementary(axisNull, other.axis)

    def __repr__(self):

        axis = "%d-Form" % len(self.axis)
        form = super().__repr__()

        return "AxialSubspace(%s, %s)" %(axis, form)

    def _initAxialSubspace(self, axis, refPt, spanVecs):

        if(not isinstance(axis, np.ndarray)): axis = np.asarray(axis)
        axis = axis.reshape((-1, EUCSPACEDIM))

        if(not isinstance(spanVecs, np.ndarray)): 
            spanVecs = np.asarray(spanVecs)
        
        self.axis = axis
        self.degree = len(self.axis)
        self.isTrans = np.all(abs(axis) <= EPSILON)
        if(self.isTrans):
            spanVecsDirected = np.asarray([[]])
        else:
            if(self.isPrinciple):
                # axis is always aligned to the basis vectors
                spanVecsDirected = axis
            else:
                spanVecsDirected = spanVecs
        
        if(self.isPrinciple):
            for vec in spanVecs:
                spanVecsDirected = Subspace.expand(spanVecsDirected, vec)
        
        super().__init__(spans=spanVecsDirected, refPt=refPt)
    
    def printInfo(self):

        print(self.__repr__())
        print("Axii: ")
        print(self.axis)

        super().printInfo(printRefPt=True, printSpan=True)

    def axisNormal(self):

        return null_space(self.axis).T

    def getSubspace(self):

        return self.copy()

    def intersect(self, other, isFreedom=True):

        assert isinstance(other, AxialSubspace),\
               "Type error: other must be a AxialSubspace"

        selfSubspace = self.getSubspace()
        otherSubspace = other.getSubspace()
        
        intersection = selfSubspace.intersect(otherSubspace)
        
        #intersection = self.intersect(other)


        if(intersection == None): return None

        if(isFreedom):
            axisCombined = self.axis
            for axis in other.axis:
                axisCombined = Subspace.expand(axisCombined, axis)
        else:
            if(self.degree > self.spanRank):
                axisCombined = other.axis
            else:
                axisCombined = intersection.spanSpace

        intAxialSubspace = AxialSubspace(axisCombined, 
                                         intersection.refPt, 
                                         intersection.spanSpace,
                                         isPrinciple=False)
        
        return intAxialSubspace

    def curScrewSpace(self, flexures):

        goalSpace = self.outputScrewSpace()

        # get intersections
        allSVs = []
        for f in flexures:
            fAS = f.axialSubspace()
            intersection = self.intersect(fAS, isFreedom=False)
            intSVs = intersection.outputScrewSpace()
            allSVs += [intSVs]
        
        if(len(allSVs) != 0):
            allSVs = np.concatenate(allSVs, axis=0)
        
        cur = np.zeros(0)
        for vec in allSVs:
            cur = Subspace.expand(cur, vec)

        return cur

    def reposition(self, point):
        
        if(not isinstance(point, np.ndarray)):
            point = np.asarray(point)
        
        if(self.isTrans): # translation is spatailly invariant
            self.refPt = point
        else: # rotation is position-dependent
            cp = self.closestPoint(point)
            
            self.refPt = cp

    def isAxisIncludedBy(self, other):

        if(other.degree == other.dim):
            # other's axis is full-ranked
            return True
        else:
            # none full rank cases

            for vec in self.axis:
                correction = Subspace.project(vec, other.axis)
                residual = vec - correction
                if(np.any(abs(residual) > EPSILON)): return False
            
            return True
    
    def isSpaceIncludedBy(self, other):
        # other has more axis and lower dimensional subspace than self
        
        # step 1: make copies to aovid aliasing
        selfDummy = copy.deepcopy(self) # avoid aliasing
        otherDummy = copy.deepcopy(other)
        
        # step 2: remove principle axis from self
        #         move refPt position
        
        if(selfDummy.isPrinciple):
            for vec in selfDummy.axis: # move everytime unwedge
                pointer = otherDummy.refPt - selfDummy.refPt
                pointerProj = Subspace.project(pointer, vec.reshape(1, -1))
                
                selfDummy.move(pointerProj)
                #selfDummy.removeSpan(pointerProj)
                selfDummy.removeSpan(vec)
        
        # step 3: remove principle axis from other
        #         do not move refPt position
        if(otherDummy.isPrinciple):
            for vec in otherDummy.axis:
                otherDummy.removeSpan(vec)

        # step 4: check if erased subspace is included by other
        isSubset = selfDummy.isSubspaceOf(otherDummy)
        
        return isSubset

    def isIncludedBy(self, other):
        
        # check if axis is included by other
        if(not self.isAxisIncludedBy(other)): 
            return False
        # at this point, self will have equal or fewer axii than other
        
        # check if subspace is included by other
        if(not self.isSpaceIncludedBy(other)): 
            return False
        
        return True

    def alignSpace(self, bases):

        # align axii
        factor = np.zeros(len(bases))

        if(not self.isPoint):
            for vec in self.spanSpace:
                dotProd = np.dot(bases, vec)
                factor += abs(dotProd)
            
            aligned = []
            for i in range(len(bases)):
                if(factor[i] > EPSILON):
                    aligned += [bases[i]]
            
            aligned = np.stack(aligned)

            self.spanSpace = aligned

    def output(self):
        
        motionType = 't' if self.isTrans else 'r'
        axis = tuple([tuple(vec) for vec in self.axis])
        spans = tuple([tuple(vec) for vec in self.spanSpace])
        try:
            refPt = tuple(self.refPt)
        except:
            refPt = (0, 0, 0)
        report = (motionType, axis, refPt, spans)

        return report
    
    def outputScrewSpace(self):

        if(self.spanRank != 0):
            anchors = self.refPt + self.spanSpace
        else:
            anchors = self.refPt.reshape((1, -1))
        
        screws = np.zeros(0)
        for axis in self.axis:
            base = screwVec(axis, self.refPt)
            screws = Subspace.expand(screws, base)
            for anchor in anchors:
                screw = screwVec(axis, anchor)
                screws = Subspace.expand(screws, screw)
                
        return screws
