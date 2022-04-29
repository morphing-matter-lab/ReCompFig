import numpy as np
from numpy.lib.twodim_base import mask_indices
from scipy.optimize.optimize import vecnorm
import sympy as sp
from scipy.linalg import null_space
from itertools import combinations

from parameters import EPSILON, EUCSPACEDIM, MAXDOFS
from messages import MSG_UNIMPLEMENTED, MSG_INVALID_GOALS, MSG_VALID_GOALS, MSG_GOALS_DONE, MSG_CHECKLIST
from utils import fixFloat
from screwVectors import lineToScrewVec, swapOperator, normScrew, screwVec, screwVecLen, decompose
from subspace import Subspace, AxialSubspace

class RotationDOF(object):
    def __init__(self, vec, refPt=None):
        self.vec = None
        self.axis = None
        self.refPt = None

        if(len(vec) == MAXDOFS):
            self.axis, mag, _, self.refPt = decompose(vec)
            self.vec = vec / mag
        else:
            self.axis = vec / np.linalg.norm(vec)
            self.refPt = refPt
            self.vec = screwVec(self.axis, self.refPt)
    
    def subspace(self):

        space = Subspace(spans=self.axis.reshape((1, -1)), refPt=self.refPt)

        return space

    def forceToMoment(self, force):
        
        axisSpace = self.subspace()
        forceSpace = force.subspace()
        
        lever = axisSpace.minDistTo(forceSpace)
        amp = lever * force.amp

        reaction = Moment(self.axis, self.refPt, amp)

        return reaction
    
    def momentToForce(self, moment, actingPt):
        
        # check axis
        axisSpace = self.subspace()
        momentSpace = moment.subspace()

        axisDist = axisSpace.minDistTo(momentSpace)
        if(axisDist < EPSILON):

            cp = axisSpace.closestPoint(actingPt)
            leverVec = actingPt - cp
            lever = np.linalg.norm(leverVec)
            leverVecUnit = leverVec / lever
            
            vec = np.cross(self.axis, leverVecUnit)
            amp = moment.amp / lever
        else:
            vec = np.zeros(EUCSPACEDIM)
            amp = 0
        
        reaction = Force(vec, actingPt, amp)

        return reaction

    def inactSpace(self, viewCenter):
        
        vec = self.vec.reshape((1, -1))
        reciprocal = DOFSystem.freeToCon(vec)
        conSpace = DOFSystem(reciprocal, isConstraint=True)
        conSpace.recenterSpaces(viewCenter)
        nullSpace = conSpace.output()[-1]
        
        return nullSpace

    def isReactive(self, actuation):

        reaction = self.reaction(actuation)

        return abs(reaction) >= EPSILON
    
    def reaction(self, actuation):
        
        if(isinstance(actuation, RotationalAct) or isinstance(actuation, Moment)):
            axisSpace = self.subspace()
            momentSpace = actuation.subspace()
            axisDist = axisSpace.minDistTo(momentSpace)
            
            if(axisDist < EPSILON):
                sign = np.dot(self.axis, actuation.axis)
                result = sign * actuation.amp
            else:
                result = 0
        
        elif(isinstance(actuation, TranslationAct) or isinstance(actuation, Force)):
            moment = self.forceToMoment(actuation)
            result = moment.amp
        
        return result

class TranslationDOF(object):
    def __init__(self, vec):

        if(len(vec) == MAXDOFS):
            self.axis, mag, _, _ = decompose(vec)
            self.vec = vec / mag
        else:
            self.axis = vec / np.linalg.norm(vec)
            self.vec = screwVec(trans=self.axis)
    
    def inactSpace(self, viewCenter):
        
        vec = self.vec.reshape((1, -1))
        
        reciprocal = DOFSystem.freeToCon(vec)
        conSpace = DOFSystem(reciprocal, isConstraint=True)
        conSpace.recenterSpaces(viewCenter)
        nullSpace = conSpace.output()[0]

        return nullSpace
    
    def isReactive(self, actuation):

        reaction = self.reaction(actuation)

        return abs(reaction) >= EPSILON
    
    def reaction(self, actuation):

        if(isinstance(actuation, RotationalAct)):
            result = 0
        elif(isinstance(actuation, TranslationAct)):
            result = np.dot(self.axis, actuation.vec)

        return result

class Force(object):
    def __init__(self, vec, actingPt, amp=None):
        
        self.vec = vec
        self.actingPt = actingPt
        self.amp = amp

    def subspace(self):

        space = Subspace(self.unitVec.reshape((1, -1)), self.actingPt)

        return space

    def toMoment(self, rotDOF):
        # axis is a screw vector

        return rotDOF.forceToMoment(self)

class Moment(object):
    def __init__(self, axis, refPt, amp=None):

        self.axis = axis
        self.refPt = refPt
        self.amp = amp

    def subspace(self):

        space = Subspace(self.axis.reshape((1, -1)), self.refPt)

        return space
    
    def toForce(self, rotDOF, actingPt):

        return rotDOF.momentToForce(self, actingPt)

class Actuation(object):
    def __init__(self, id):
        self.id = id
        self.amp = None
        self.actingStage = None
    
    def subspace(self):

        assert False, MSG_UNIMPLEMENTED

        return 42

    def attachTo(self, stage):

        self.actingStage = stage
        stage.addActuation(self)

    @staticmethod
    def model(info):

        objType = info[0]
        for actType in Actuation.__subclasses__():
            if objType == actType.name:
                modeled = actType(info)

        return modeled

class TranslationAct(Actuation):
    name = "ActTrans"

    def __init__(self, info):
        super().__init__(info[1])

        self.actingPt = None
        self.unitVec = None
        self.amp = None
        self.vec = None

        startPt = np.asarray(info[2]).astype(np.float)
        endPt = np.asarray(info[3]).astype(np.float)

        self.actingPt = startPt
        vec = endPt - startPt
        vecNorm = np.linalg.norm(vec)

        if(vecNorm >= EPSILON):
            self.unitVec = vec / vecNorm
        else:
            self.unitVec = vec
        
        self.amp = info[4]
        self.vec = self.unitVec * self.amp

    def subspace(self):

        space = Subspace(self.unitVec.reshape((1, -1)), self.actingPt)

        return space

    def toMoment(self, rotDOF):
        # axis is a screw vector

        return rotDOF.forceToMoment(self)

class RotationalAct(Actuation):
    name = "ActRot"

    def __init__(self, info):
        super().__init__(info[1])

        self.refPt = None
        self.axis = None
        self.amp = None

        startPt = np.asarray(info[2]).astype(np.float)
        endPt = np.asarray(info[3]).astype(np.float)

        self.refPt = startPt
        self.amp = info[4]

        axis = endPt - startPt
        self.axis = axis / np.linalg.norm(axis)

    def subspace(self):

        space = Subspace(self.axis.reshape((1, -1)), self.refPt)

        return space
    
    def toForce(self, rotDOF, actingPt):

        return rotDOF.momentToForce(self, actingPt)

class DOFSystem(object):
    def __init__(self, twistVecs, isConstraint=False):
        
        self.model = None

        self.tVecs = twistVecs
        self.deg = None
        self.nullRank = None
        self.nullSpace = None
        self.bases = None

        self.wrenches = None
        self.coords = np.asarray(sp.symbols("x y z"))

        self.permTrans = []
        self.permRots = []
        self.freedoms = []
        self.motionSubspaces = []

        self.isConstraint = isConstraint

        self._resolveSystem()

    def _resolveSystem(self):

        self.tVecs = DOFSystem.cleanVecs(self.tVecs)
        self._checkVecs()
        
        self._getWrenchParameters()
        self._findTranslations()
        self._findBases()
        self._findRotations()
        self._updateMotionSpaces()

    def _checkVecs(self):

        norms = np.linalg.norm(self.tVecs, axis=-1)
        
        valid = []
        for i in range(len(self.tVecs)):
            if(norms[i] > EPSILON):
                normVec = self.tVecs[i] / screwVecLen(self.tVecs[i])
                valid += [normVec]
        
        self.tVecs = np.asarray(valid)
        self.deg = self.tVecs.shape[0]

    def _findBases(self):
        
        if(self.deg == 0):
            self.bases = np.identity(EUCSPACEDIM)
            return

        hasTrans = len(self.permTrans) > 0
        # find rotation vector subsapce
        axii = self.tVecs[:,:EUCSPACEDIM]
        axisSubspace = Subspace(spans=axii, parameters=self.coords)

        
        if(hasTrans):
            # find translation subspace
            trans = self.permTrans[0][2].output(self.tVecs, mode='t')
            tranVecs = DOFSystem.cleanVecs(np.asarray(trans["spanVecs"]))
            transSubspace = Subspace(spans=tranVecs, parameters=self.coords)
            # find "unique" vector in the two subspaces and their intersection
            spanAxisInt = transSubspace.intersect(axisSubspace)
            spacesToCheck = [transSubspace, axisSubspace, spanAxisInt]
        else:
            # rotation only
            spacesToCheck = [axisSubspace]
            
        principleVecs = []
        for subspace in spacesToCheck:
            if(subspace.normalRank == 1):
                principleVecs += [subspace.normalSpace]
            elif(subspace.spanRank == 1):
                principleVecs += [subspace.spanSpace]

        if(len(principleVecs) == 0): # rotation in all directions or none at all
            self.bases = np.eye(EUCSPACEDIM)
        else:
            principleVecs = np.concatenate(principleVecs, axis=0)
            principleSubspace = Subspace(spans=principleVecs, 
                                         refPt=[0,0,0], 
                                         parameters=self.coords)
            self.bases = principleSubspace.getBases()
        
        # check if axis aligned
        aligned = True
        for axis in self.bases:
            axis = axis / np.linalg.norm(axis)
            hots = np.where(np.abs(axis) > EPSILON)
            if(len(hots[0]) != 1):
                aligned = False
        
        if(aligned):
            self.bases = np.identity(EUCSPACEDIM)
        
        """
        # special case, if bases are somewhat axes aligned, then use world-space axes
        aligned = False
        for axis in self.bases:
            axis = axis / np.linalg.norm(axis)
            hots = np.where(np.abs(axis) > EPSILON)
            if(len(hots[0]) == 1):
                aligned = True
        
        if(aligned):
            self.bases = np.identity(EUCSPACEDIM)
        """
         
    def _getWrenchParameters(self):

        if(self.deg == 0):
            self.wrenches = []
            return

        wrenchSymbols = " ".join(["w%d" % (i + 1) for i in range(self.deg)])
        wrenches = sp.symbols(wrenchSymbols)
        if(self.deg == 1): wrenches = [wrenches]

        self.wrenches = wrenches

    def _findTranslations(self):

        if(self.deg == 0):
            self.permTrans = []
            return

        axisTarget = np.zeros(EUCSPACEDIM)
        vecs = self.tVecs[:, :EUCSPACEDIM]
        permissibleTrans = []
        system = np.concatenate((vecs.T, axisTarget.reshape(-1, 1) * -1), axis=-1)
        solution = self._findSolutionSpace(system, self.wrenches)
        if(solution != None and\
           not (solution.isPoint and\
                abs(np.linalg.norm(solution.refPt)) <= EPSILON)):
            permissibleTrans += [('t', axisTarget, solution)]
            
        self.permTrans = permissibleTrans
        
    def _findRotations(self):
        
        if(self.deg == 0):
            self.permTrans = []
            return
        vecs = self.tVecs[:, :EUCSPACEDIM]
        #vecs[abs(vecs) < 1e-4] = 0
        
        trans = self.tVecs[:, EUCSPACEDIM:]
        permissibleRots = []
        for basis in self.bases:
            #basis[abs(basis) < 1e-4] = 0
            
            # pure rotation constraint
            zeroTransLeft = trans.T * basis.reshape(-1, 1)
            zeroTrans = np.sum(zeroTransLeft, axis=0)
            zeroTrans = np.concatenate((zeroTrans, np.zeros(1))).reshape(1, -1)
            # construct system of equations
            system = np.concatenate((vecs.T, basis.reshape(-1, 1) * -1), axis=-1)
            system = np.concatenate((system, zeroTrans), axis=0)
            
            solutions = self._findSolutionSpace(system, self.wrenches)

            if(solutions != None):
                permissibleRots += [('r', basis, solutions)]
                
           
        self.permRots = permissibleRots

    def _findSolutionSpace(self, system, wrenches):

        dim = len(wrenches)
        wrenches = np.asarray(wrenches)

        subspaces = []
        for i in range(len(system)):
            row = system[i]
            vec = row[:-1] # excluding constants
            offset = row[-1] * -1
            vecLen = np.linalg.norm(vec)
            if(vecLen <= EPSILON and abs(offset) > EPSILON): 
                return None # shortcut: impossible solution
            if(vecLen <= EPSILON and abs(offset) <= EPSILON): 
                continue # shortcut: immediately evaluated to True

            vecLenFixed = vecLen if abs(vecLen) >= EPSILON else 1 # avoid div-0
            vecUnitized = vec / vecLenFixed # normal space is support vector
            supportVecLen = offset / vecLenFixed
            refPt = vecUnitized * supportVecLen # support vector as ref point

            subspace = Subspace(normals=vecUnitized, refPt=refPt, 
                                parameters=wrenches)
            subspaces += [subspace]
        
        subspaces = Subspace.cullDuplicates(subspaces)
        
        if(len(subspaces) == 0):
            solutions = Subspace(refPt=np.zeros(dim), spans=np.identity(dim),
                                 parameters=wrenches)
        else:
            solutions = Subspace.solve(subspaces)
            
        return solutions

    def _updateMotionSpaces(self):

        self.freedoms = self.permRots + self.permTrans
        
        self.rotDeg = len(self.permRots)
        self.transDeg = self.deg - self.rotDeg
        self._findMotionSpaces()
        self.alignSpaces()

    def _findMotionSpaces(self):

        trans, rots = [], []
        
        for freedom in self.freedoms:
            
            # (type, axis, dict["refPt", "spanVecs"])
            motionType, axis, subspace = freedom
            subspace = subspace.output(self.tVecs, mode=motionType)
            motionSubspace = AxialSubspace(axis, subspace["refPt"], subspace["spanVecs"])

            if(motionType == 'r'): rots += [motionSubspace]
            elif(motionType == 't'): trans += [motionSubspace]
            else: assert False, MSG_UNIMPLEMENTED
            
        rots = self._findAllRotSpace(rots)
        self.motionSubspaces = rots + trans

    def _findAllRotSpace(self, rots):
        
        basesRotsCount = len(rots)
        newRots = []
        # 2-rotation DOF
        if(basesRotsCount >= 2):
            combs = combinations(rots, 2)
            for r1, r2 in combs:
                
                intersection = r1.intersect(r2)
                
                if(intersection != None):
                    newRots += [intersection]
        
        # 3-rotation DOF
        if(basesRotsCount >= 3):
            intFound = True
            intersection = rots[0]
            for rot in rots[1:]:
                intersection = intersection.intersect(rot)
                if(intersection == None): 
                    intFound = False
                    break
            if(intFound): newRots += [intersection]
        
        rots = rots + newRots
        uniques = []
        for i in range(len(rots)):
            isUnique = True
            for j in range(len(rots)):
                if(i == j): continue
                r1, r2 = rots[i], rots[j]
                isIncluded = r1.isIncludedBy(r2)
                if(isIncluded):
                    isUnique = False
                    break
            if(isUnique):
                uniques += [i]
        rots = [rots[i] for i in uniques]
        
        return rots

    def solveTwist(self, motion):

        colVec = motion.reshape((-1, 1))
        system = np.concatenate((self.tVecs.T, colVec * -1), axis=-1)
        solution = self._findSolutionSpace(system, self.wrenches)
        
        return solution

    def solveTrans(self, axis):

        rotTarget = np.zeros(EUCSPACEDIM).reshape((-1, 1))
        axis = axis.reshape((-1, 1))
        target = np.concatenate((rotTarget, axis), axis = 0)
        system = np.concatenate((self.tVecs.T, target * -1), axis=-1)
        solution = self._findSolutionSpace(system, self.wrenches)

        return solution
    
    def solveRot(self, axis):

        vecs = self.tVecs[:, :EUCSPACEDIM]
        trans = self.tVecs[:, EUCSPACEDIM:]

        # pure rotation constraint
        zeroTransLeft = trans.T * axis.reshape(-1, 1)
        zeroTrans = np.sum(zeroTransLeft, axis=0)
        zeroTrans = np.concatenate((zeroTrans, np.zeros(1))).reshape(1, -1)
        # construct system of equations
        system = np.concatenate((vecs.T, axis.reshape(-1, 1) * -1), axis=-1)
        system = np.concatenate((system, zeroTrans), axis=0)

        solution = self._findSolutionSpace(system, self.wrenches)

        return solution

    def printInfo(self):

        axii = np.take(self.tVecs, range(EUCSPACEDIM), axis=-1)
        pts = np.take(self.tVecs, range(EUCSPACEDIM, MAXDOFS), axis=-1)
        print("")
        print("axii")
        print(axii)
        print("")

        print("points")
        print(pts)
        print("")
        print("system basis")
        print(self.bases)

    def findSolution(self, motion):
        # motion is a screw vector
        if(isinstance(motion, tuple) and len(motion) == 3):
            targetMotion = lineToScrewVec(motion)
        else:
            targetMotion = motion
        
        targetMotion = normScrew(targetMotion)
        solution = self.solveTwist(targetMotion)
        if(solution == None): return None
        else:
            sol = solution.outputRefPtSpan()
            return sol

    def isMotionPermissible(self, motions):
        
        result = []
        
        for motion in motions:
            targetMotion = lineToScrewVec(motion)
            
            vecNorm = np.linalg.norm(targetMotion)
            vec = targetMotion[:EUCSPACEDIM]

            if(vecNorm <= EPSILON): # no motion
                result += [True]
                continue

            solution = self.solveTwist(targetMotion)
            
            if(solution == None): result += [False]
            else: result += [True]
        
        return result

    def freedom(self, printResults= False):

        if(printResults):
            print()
            print("wrench parameters:")
            print(self.wrenches)

            print("permitted rotations:")
            if(len(self.permRots) == 0):
                print("None")
            else:
                for item in self.permRots:
                    print(item)
            
            print("permitted translations:")
            if(len(self.permTrans) == 0):
                print("None")
            else:
                for item in self.permTrans:
                    print(item)
            print()

        return self.freedoms

    def recenterSpaces(self, viewPt):
        
        if(viewPt != None):
            for motionSubspace in self.motionSubspaces:
                motionSubspace.reposition(viewPt)

    def alignSpaces(self):

        for motionSubspace in self.motionSubspaces:
            motionSubspace.alignSpace(self.bases)

    def output(self, axialSubspaceMode=False):

        spaces = []

        for motionSubspace in self.motionSubspaces:
            
            if(self.isConstraint and motionSubspace.isTrans): continue
            
            if(axialSubspaceMode):
                result = motionSubspace
            else:
                result = motionSubspace.output()
            
            spaces += [result]
        
        return spaces

    def outputBases(self):

        bases = tuple([tuple(vec) for vec in self.bases])
        
        return bases

    def constraintSys(self):

        constraintVecs = DOFSystem.freeToCon(self.tVecs)
        constraintSys = DOFSystem(constraintVecs, isConstraint=True)

        return constraintSys

    def isAxial(self):

        axial = self.tVecs[:,:EUCSPACEDIM]
        isValid = np.linalg.norm(axial) > EPSILON

        return isValid

    def isEmpty(self):

        size = self.tVecs.size
        
        return size == 0

    def constraintDelta(self, tarMotions, viewCenter):

        # figure out target motion's required constraint space
        freedom = np.stack([lineToScrewVec(tm) for tm in tarMotions])
        consNeeded = DOFSystem.freeToCon(freedom)
        #check if target motion is achievable
        valid, goalsInfo = DOFSystem.checkGoalValidity(consNeeded)
        
        if(not valid): return None, goalsInfo

        # check overconstraint here
        isOverCon, overConsInfo = self.model.checkOverConstraint()
        
        # generate header message
        header = goalsInfo + '\n' + overConsInfo + '\n'

        if(not isOverCon):
            # at this point, the target DOF is achievable with the current setup
            spaces, msg, isDone = self.checkConsSubspaces(consNeeded, viewCenter)
            if(isDone):
                header = MSG_GOALS_DONE + '\n' + overConsInfo + '\n'
            msg = [header + m for m in msg]
        else:
            spaces, msg = [], []
        
        return spaces, msg

    def checkConsSubspaces(self, consNeeded, viewCenter):
        
        conSys = DOFSystem(consNeeded, True)
        
        conSys.model = self.model
        conSys.recenterSpaces(viewCenter)
        flexInGroup = conSys.model.getPhJointCons()
        subspaces = conSys.output(axialSubspaceMode=True)
        deltaSubspaces = []
        checks = []

        isCompleted = True
        for subspace in subspaces:
            
            goal = subspace.outputScrewSpace()
            
            tarAxisDeg = subspace.degree
            tarSpanDeg = subspace.spanRank
            minWires = len(null_space(null_space(goal).T).T)
            # find flexure elements associated with subspace
            associatedCons = []
            # association: intersects at more than a point
            for f in flexInGroup:
                isAssociated = f.isAssociated(subspace)
                
                if(isAssociated):
                    associatedCons += [f]
            
            cur = subspace.curScrewSpace(associatedCons)
            deltaSpace = DOFSystem.deltaSpace(goal, cur)
            
            deltaSubspaces += [subspace.output()]
            
            curSys = DOFSystem(cur)
            curSys.recenterSpaces(viewCenter)
            deltaSys = DOFSystem(deltaSpace)
            deltaSys.recenterSpaces(viewCenter)
            
            if(len(deltaSpace) != 0):
                isCompleted = False

            curSysSubspace = curSys.output(axialSubspaceMode=True)
            if(len(curSysSubspace) == 0):
                curAxisDeg = 0
                curSpanDeg = 0
            else:
                curAxisDeg = max([s.degree for s in curSysSubspace])
                curSpanDeg = max([s.spanRank for s in curSysSubspace])
            curWires = len(cur)
            
            msg = MSG_CHECKLIST(tarAxisDeg, tarSpanDeg, minWires,\
                                curAxisDeg, curSpanDeg, curWires)

            checks += [msg]
        
        return deltaSubspaces, checks, isCompleted

    def orthoDOF(self):
        
        vecs = []
        # returns more than basics
        for f in self.motionSubspaces:
            
            if(f.isTrans):
                for vec in f.spanSpace:
                    vecs += [screwVec(trans=vec)]
            else:
                for vec in f.axis:
                    vecs += [screwVec(axis=vec, ref=f.refPt)]
        
        if(len(vecs) > 0):
            vecs = np.stack(vecs)
        else:
            vecs = np.zeros(0)
            
        return vecs

    def basesRots(self, center):
        
        return np.stack([screwVec(vec, center) for vec in self.bases])
  
    def viewInactSpace(self, viewCenter):

        spaces = []
        for freedom in self.freedoms:
            
            motionType, axis, subspace = freedom
            subspace = subspace.output(self.tVecs, mode=motionType)
            if(motionType == 'r'):
                DOF = RotationDOF(axis, subspace["refPt"])
                space = DOF.inactSpace(viewCenter)
                spaces += [space]
            elif(motionType == 't'):
                for axis in subspace["spanVecs"]:
                    DOF = TranslationDOF(axis)
                    space = DOF.inactSpace(viewCenter)
                    spaces += [space]
        
        return spaces

    def checkActuation(self, act, threshold=EPSILON):

        result = []
        for freedom in self.freedoms:
            
            motionType, axis, subspace = freedom
            subspace = subspace.output(self.tVecs, mode=motionType)
            if(motionType == 'r'):
                DOF = RotationDOF(axis, subspace["refPt"])
                reaction = DOF.reaction(act)
                actuated = abs(reaction) > threshold
                result += [actuated]
            elif(motionType == 't'):
                for axis in subspace["spanVecs"]:
                    DOF = TranslationDOF(axis)
                    reaction = DOF.reaction(act)
                    actuated = abs(reaction) > threshold
                    result += [actuated]
                
        return result

    def findFreeRotations(self):
        
        freeRots = []
        vecs = self.orthoDOF()
        for vec in vecs:
            if(np.linalg.norm(vec[:EUCSPACEDIM]) > EPSILON):
                freeRots += [vec]
        
        freeRots = np.stack(freeRots)

        return freeRots

    def indiDOFs(self):
        
        result = {}
        axesName = ['x', 'y', 'z']
        
        result = {}
        for freedom in self.freedoms:

            motionType, axis, subspace = freedom
            # rotation
            if(motionType == 'r'):
                # find id
                axesDot = np.dot(axis, self.bases)
                id = np.where(axesDot != 0)[0][0]
                others = np.concatenate([self.bases[:id], self.bases[id + 1:]])

                rotVec = [screwVec(axis, np.zeros(EUCSPACEDIM))]
                transVec = [screwVec(trans=vec) for vec in others]
                rotFreedom = np.stack(rotVec + transVec)

                tag = 'r' + axesName[id]
                result[tag] = rotFreedom

            elif(motionType == 't'):
                subspace = subspace.output(self.tVecs, mode=motionType)
                vecs = subspace["spanVecs"]
                for vec in vecs:
                    # find id
                    axesDot = np.dot(vec, self.bases)
                    id = np.where(axesDot != 0)[0][0]
                    
                    transFreedom = screwVec(trans=vec).reshape(1, -1)
                    
                    tag = 't' + axesName[id]
                    result[tag] = transFreedom
            
        return result

    @staticmethod
    def normVecs(vecs):


        length = np.linalg.norm(vecs[:,:EUCSPACEDIM], axis=1).reshape((-1, 1))
        # replace 0 with 1 to avoid divide by zero errors
        length[np.isclose(length, np.zeros(length.shape))] = 1

        vecs = vecs / length

        return vecs

    @staticmethod
    def cleanVecs(vecs):
        
        new = []
        for vec in vecs:
            vecLen = screwVecLen(vec)
            v = normScrew(vec)
            v = fixFloat(v)
            v = v * vecLen
            new += [v]

        if(len(new) == 0):
            return vecs
        result = np.stack(new)
        #vec = DOFSystem.normVecs(vec)
        #vec = fixFloat(vec)
        #vec = DOFSystem.normVecs(vec)
        
        return result

    @staticmethod
    def conToFree(con):
        
        I = swapOperator()
        
        right = np.matmul(I, con.T)
        #freedom = fixFloat(null_space(right.T).T)
        freedom = null_space(right.T).T
        
        return freedom

    @staticmethod  
    def freeToCon(free):

        I = swapOperator()
        left = np.matmul(free, I)
        
        #constraint = fixFloat(null_space(left).T)
        constraint = null_space(left).T

        return constraint

    @staticmethod
    def test(tarMotions):
        
        goals = [lineToScrewVec(motion) for motion in tarMotions]
        goals = np.asarray(goals)
        goalsConSpace = DOFSystem.freeToCon(goals)
        print(goals)
        print(goalsConSpace)

    @staticmethod
    def checkGoalValidity(consSpace):

        isEmptyCons = consSpace.size == 0
        isAxisless = np.linalg.norm(consSpace[:,:EUCSPACEDIM]) < EPSILON
        
        if(isEmptyCons or isAxisless): # unachievable
            return False, MSG_INVALID_GOALS
        else:
            return True, MSG_VALID_GOALS
    
    @staticmethod
    def spaceDiff(first, second):

        # generate header message
        header = goalsInfo + '\n' + overConsInfo + '\n'

        if(not isOverCon):
            # at this point, the target DOF is achievable with the current setup
            spaces, msg, isDone = self.checkConsSubspaces(consNeeded, viewCenter)
            if(isDone):
                header = MSG_GOALS_DONE + '\n' + overConsInfo + '\n'
            msg = [header + m for m in msg]
        else:
            spaces, msg = [], []
        
        return spaces, msg


    @staticmethod
    def deltaSpace(space1, space2):
        
        space = np.copy(space1)
        for vec in space2:
            space = Subspace.remove(space, vec)
        
        return space

if __name__ == "__main__":
    _1 = [0, 0, 1, 0, 0, 0]
    _2 = [0, 1, 0, 0, 0, 0]
    _3 = [1, 0, 0, 0, 0, 0]
    _4 = [0, 0, 0, 0, 0, 1]
    _5 = [0, 0, 0, 0, 1, 0]
    _6 = [0, 0, 0, 1, 0, 0]
    f = [_1, _2]
    #c = [_1, _2, _4, _5, _6]
    #c = np.asarray(c)
    c = DOFSystem.freeToCon(f)
    sys = DOFSystem(c, True)
    o = sys.output()
    print(o)