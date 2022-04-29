from decimal import ConversionSyntax
from types import SimpleNamespace
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import null_space
from scipy.optimize._trustregion_constr.qp_subproblem import modified_dogleg
from scipy.optimize.nonlin import NoConvergence
from itertools import combinations

from scipy.sparse import construct

from parameters import EUCSPACEDIM, MAXDOFS, EPSILON, PRINT_NUM, ROD_MIN_ASPECT_RATIO, ROD_MAX_ASPECT_RATIO
from screwVectors import Frame, lineToScrewVec, screwTransform, screwVec, getAxis
from utils import fixFloat, str2num
from DOFsystem import Actuation, DOFSystem
from subspace import Subspace
from messages import *
from graph import Graph
from joints import Flexure, Imaginary, PlaceHolder
from stages import Stage
from material import Material

from input import*

class CompMech(object):
    def __init__(self, inputModel):

        self.nodeCount = None
        self.edgeCount = None
        self.A = None # adjacency matrix
        self.C = None # incidency matrix
        self.Q = None # fewest closed loop matrix

        self.flexures = None
        self.joints = None
        self.stages = None
        self.materials = None
        self.goals = None
        self.actuation = None

        self.jfts = None
        self.FT = None # freedom topology
        self.X = None # system DOF

        self.ground = None
        self.target = None
        self.P = None # path vector

        self.frame = None
        
        self.parameters = None

        self._modelSystem(inputModel)

    def _modelSystem(self, model):
        
        model = self._preproc(model)
        self._makeModel(model)
        self._makeGraph()
        self._makeJfts()
        self._makeX()

        #self.solveDeformation()

    def _preproc(self, model):
        
        def str2num(inObj):

            if(isinstance(inObj, tuple) or isinstance(inObj, list)):
                container = []
                for obj in inObj:
                    container += [str2num(obj)]
                
                if(isinstance(inObj, tuple)): # convert to same type as input
                    container = tuple(container)
                return container
            else:
                if(isinstance(inObj, str)):
                    try:
                        return float(inObj)
                    except:
                        return inObj
                else:
                    return inObj

        converted = str2num(model)

        return converted

    def _makeGraph(self):
        
        # make dictionary form of the graph
        dictForm = {}
        dictForm["nodeCount"] = len(self.stages)
        dictForm["edgeCount"] = len(self.flexures)
        
        adjDict = {}
        for flexure in self.flexures:
            id, adj = flexure.id, flexure.adj
            adjDict[id] = adj
        dictForm["adjacency"] = adjDict

        self.graph = Graph(dictForm)
        
    def _makeModel(self, model):

        # frame
        self.frame = Frame.spatial()
        # hardcoded word coord

        # materials
        self.materials = []
        for i in range(len(model[1])):
            matInfo = model[1][i]
            material = Material.create(matInfo)
            self.materials += [material]
        
        # flexures
        self.flexures = []
        for i in range(len(model[2])):
            flexureDesign = model[2][i]
            matId = model[2][i][3]
            material = self.materials[matId]
            flexure = Flexure.model(flexureDesign, material)
            self.flexures += [flexure]
        self._addImaginaryFlexures()

        # TODO: need to model joints
        
        # stages
        self.stages = []
        for i in range(len(model[4])):
            stageDesign = model[4][i]
            matId = model[4][i][2]
            material = self.materials[matId]
            stage = Stage.model(stageDesign, material)
            self.stages += [stage]

        # goals
        self.goals = []
        for i in range(len(model[5])):
            goal = model[5][i]
            self.goals += [goal]
        
        # chain joints and stages
        for joint in self.flexures:
            baseId, targId = joint.adj
            base, targ = self.stages[baseId], self.stages[targId]
            joint.chain(base, targ)

        self.actuation = []
        for i in range(len(model[6])):
            info = model[6][i]
            actingStage = self.stages[info[5]]
            actuation = Actuation.model(info)
            actuation.attachTo(actingStage)
            self.actuation += [actuation]

        """
        for joint in self.flexures:
            joint.subDivide()
        """

    def _addImaginaryFlexures(self):
        adjDict = {}
        for flexure in self.flexures:
            adjDict[flexure.adj] = adjDict.get(flexure.adj, []) + [flexure]
        
        for adj in adjDict:
            flexures = adjDict[adj]
            if(len(flexures) == 1):
                id = len(self.flexures)
                material = flexures[0].material
                info = (Imaginary.name, id, adj)
                flexure = Flexure.model(info, material)
                self.flexures += [flexure]

    def _makeJfts(self):
        
        jfts = []
        for f in self.flexures:
            freedom = f.freedom
            jfts += [freedom.T]
        
        self.jfts = jfts

    def _makeX(self):

        Q = self.graph.Q
        FTRows = Q.shape[0] * MAXDOFS
        FTCols = sum([jft.shape[1] for jft in self.jfts])
        
        FT = np.zeros((FTRows, FTCols))
        for row in range(Q.shape[0]):
            DOFSum = 0
            for col in range(Q.shape[1]):
                jft = self.jfts[col]
                jDOFs = jft.shape[1]
                factor = Q[row, col]
                factorized = factor * jft
                FT[row * MAXDOFS : (row + 1) * MAXDOFS, DOFSum : DOFSum + jDOFs] = factorized
                DOFSum += jDOFs

        # find null sapce of freedom topology
        FT[FT == np.NZERO] = 0
        
        if(FT.size != 0):
            #X = null_space(FT)
            
            X = Subspace(spans=FT).normalSpace.T
            
        else:
            X = np.zeros((0, 0))
        
        self.X = X

    def _makeP(self, ground, target):

        pathVec = self.graph.makeP(ground, target)

        P = np.zeros((0, 0)).reshape((MAXDOFS, -1))
        for i in range(len(pathVec)):
            factor = pathVec[i]
            jft = self.jfts[i]
            factorized = factor * jft
            P = np.append(P, factorized, axis=-1)
        
        return P

    def _solveDOF(self, P):
        
        if(self.edgeCount == 1):
            twistVecs = self.flexures[0].freedom
        else:
            twistVecs = np.matmul(P, self.X).T
            
        DOFsys = DOFSystem(twistVecs)
        DOFsys.model = self
        
        return DOFsys

    def findFreedom(self, ground, target):
        
        P = self._makeP(ground, target)
        
        DOFsys = self._solveDOF(P)

        return DOFsys

    def findMotion(self, ground, target, actuation):

        if(ground == target): return np.zeros(MAXDOFS)

        P = self._makeP(ground, target)
        
        sysAct = np.matmul(self.X, actuation)
        targetAct = fixFloat(np.matmul(P, sysAct))

        return targetAct

    def getPhJointCons(self):

        # find place holder
        phAdj = None
        for f in self.flexures:
            if(isinstance(f, PlaceHolder)):
                phAdj = f.adj
                break
        
        if(phAdj == None):
            return None

        # find flexures that belongs to the place holder
        flexInGroup = []
        for f in self.flexures:
            if(f.isReal and f.adj == phAdj):
                flexInGroup += [f]

        return flexInGroup

    def checkOverConstraint(self):
        # find goal constraint space
        freedom = np.stack([lineToScrewVec(tm) for tm in self.goals])
        consNeeded = DOFSystem.freeToCon(freedom)
        
        consPerp = Subspace(spans=consNeeded)
        consVecs = consPerp.spanSpace
        
        flexInGroup = self.getPhJointCons()
        if(flexInGroup == None): return MSG_NO_PLACEHOLDER
        
        # check if each component is allowed
        isValid = True
        invalidFlexures = []
        for f in flexInGroup:
            isAllowed = f.isAllowed(consVecs)
            if(not isAllowed):
                invalidFlexures += [f]
                isValid = False
        
        # generate message
        msg = MSG_OVER_CONS(invalidFlexures)
        
        return not isValid, msg

    def checkActuations(self, ground, target=None):

        result = []
        # check specific stage
        if(target != None):
            temp = []
            stageDOF = self.findFreedom(ground, target)
            
            for act in self.stages[target].actuations:
                checked = stageDOF.checkActuation(act)
                temp += [checked]
            
            for j in range(len(stageDOF.freedoms)):
                actuated = False
                for i in range(len(temp)):
                    actuated = actuated or temp[i][j]
                result += [actuated]
        else:
            # check all
            for act in self.actuation:
                tarStageId = act.actingStage.id
                stageDOF = self.findFreedom(ground, tarStageId)
                checked = stageDOF.checkActuation(act)
                result += [checked]
        
        return result

    def solveDeformation(self, inputWrenches=[], inputDeformation=[]):


        self._setActingWrenches(inputWrenches)
        self._setTargetDeformation(inputDeformation)

        self._gatherAllVariables()
        goalFunc = self._goalFunc()
        constraints = self._kinematicConstraints()
        constraints += [{'type': 'eq', 'fun': lambda x: goalFunc(x)}]

        initGuess = np.random.normal(np.zeros(len(self.parameters)))
        #print(np.matmul(self.flexures[0].compliance(), self.stages[1].actingWrench))
        initGuess *= 0
        sol = minimize(goalFunc, initGuess, constraints=constraints, options={"maxiter": 1000})
        
        for cons in constraints:
            print('\t', cons['type'], cons['fun'](sol.x))
        #print(sol)"""
        """
        initGuess = np.random.normal(np.zeros(len(self.parameters)))
        initGuess = np.zeros(len(self.parameters))
        print(goalFunc(initGuess))
        #print(self.parameters)
        sol = minimize(goalFunc, initGuess, constraints=constraints)#, options={"maxiter": 1000})
        print(sol)
        """
        """
        for i in range(50):
            initGuess = np.random.normal(np.zeros(len(self.parameters)))

            sol = minimize(goalFunc, initGuess, constraints=constraints)#, options={"maxiter": 1000})

            c = [cons['fun'](sol.x) for cons in constraints]
            print(sol.message)
            print(goalFunc(sol.x), c)

            self._setVariables(sol.x)
            print(decompose(self.stages[1].twistVars))
            #print(i)
        """
        pass

    def outputStageFrames(self):

        result = []
        for stage in self.stages:
            result += [stage.frame.output()]
        
        return tuple(result)

    def tarSolveFastSim(self, ground, target):

        DOFsys = self.findFreedom(ground, target)
        actuation = DOFsys.findSolution(self.goals[0])
        
        if(actuation is None):
            return None
        
        result = tuple(self.freedomToFreeScrew(actuation, ground, isActuation=True))

        return result

    def freeFastSim(self, ground, target):

        result = []
        DOFsys = self.findFreedom(ground, target)
        
        # rotations
        for f in DOFsys.freedoms:
            if(f[0] == 'r'):
                thisAct = self.freedomToFreeScrew(f, ground)
            elif(f[0] == 't'):
                continue
            result += thisAct
        
        # translations
        for axis in DOFsys.outputBases():
            trans = screwVec(trans=axis)
            act = DOFsys.findSolution(trans)
            if(not isinstance(act, np.ndarray)):
                continue

            thisAct = self.freedomToFreeScrew(act, ground, isActuation=True)
            result += thisAct
        
        result = tuple(result)
        
        return result

    def freedomToFreeScrew(self, freedom, ground, isActuation=False):

        if(isActuation):
            actVecs = freedom
        else:
            if(freedom[0] == 'r'):
                actVecs = freedom[2].outputRefPtSpan()
            elif(freedom[0] == 't'):
                actVecs = freedom[2].outputRefPtSpan(spanOnly=True)
            else:
                assert False, MSG_UNIMPLEMENTED
            
        result = []
        for stage in self.stages:
            perStage = []
            for vec in actVecs:
                transScrew = self.findMotion(ground, stage.id, vec)
                perStage += [tuple(transScrew)]
            result += [tuple(perStage)]
        
        result = [tuple(result)]
        
        return result

    def _setActingWrenches(self, wrenches):
        # for each stage in system, set an acting wrench
        s1Tar = np.asarray([-100, 0, 0, 0, 0, 0])
        s1Frame = self.stages[1].frame
        Ad = s1Frame.AdToWorld()
        s1Tar = np.matmul(Ad, s1Tar)

        self.stages[1].setWrench(s1Tar)
        pass 

    def _setTargetDeformation(self, deformations):
        # for each stage in system, set an prescribed deformation
        
        for stage in self.stages:
            stage.setFixedEnd()
        pass

    def _strainEnergy(self):

        # calculate energy
        results = []
        for i in range(len(self.flexures)):
            joint = self.flexures[i]
            for s in joint.segments:
                segStrainEnergy = s.strainEnergy()
                results += [segStrainEnergy]
        
        energy = sum(results)

        return energy
    
    def _negativeWork(self):

        results = []
        for i in range(len(self.stages)):
            stage = self.stages[i]
            negWork = stage.negWork()
            results += [negWork]
        
        energy = sum(results)
        return energy

    def _kinematicConstraints(self):

        # constraint equation builder
        def _eqMakerTwistChain(joint):

            def newFunc(x, printInfo=PRINT_NUM):
                self._setVariables(x)
                constraintVal = joint.kinematicConstraints()
                if(printInfo):
                    print("\t%.3f" % constraintVal)
                return constraintVal
            
            return newFunc

        def _eqMakerStageWork(stage):
            
            def newFunc(x, printInfo=PRINT_NUM):
                self._setVariables(x)
                work = -1 * stage.negWork() - EPSILON
                if(printInfo):
                    print("\t%.3f" % work)
                return work

            return newFunc

        def _eqZeroGrad():

            def newFunc(x, printInfo=PRINT_NUM):
                vals = np.zeros(len(x))
                for i in range(len(x)):
                    modifier = np.zeros(len(x))
                    modifier[i] += EPSILON
                    xMod = x + modifier

                    val0 = self.goalFunc(x)
                    val1 = self.goalFunc(xMod)
                    vals[i] = (val1 - val0) / EPSILON
                
                return np.linalg.norm(vals)
            
            return newFunc
                
        # get all constraints for joints
        constraints = []
        for i in range(len(self.flexures)):
            joint = self.flexures[i]
            consFunc = _eqMakerTwistChain(joint)
            cons = {'type': 'eq', 'fun': consFunc}
            constraints += [cons]
        
        consFunc = _eqZeroGrad()
        cons = {'type': 'eq', 'fun': consFunc}
        constraints += [cons]
        
        self.constraints = constraints
        return constraints

    def _gatherAllVariables(self):

        parameters = []
        for obj in self.flexures + self.stages:
            parameters += obj.getPrms()
        
        self.parameters = parameters

    def _goalFunc(self):

        def gF(x, printInfo=PRINT_NUM):
            # set variables
            self._setVariables(x)
            
            # calculate energy
            strain = self._strainEnergy()
            negWork = self._negativeWork()
            value = strain + negWork

            if(printInfo): 
                print(strain, negWork, value)

            return value
        
        self.goalFunc = gF
        return gF

    def _setVariables(self, vals, printMap=False):

        # create a map of variables
        valMap = {}
        for name, val in zip(self.parameters, vals):
            valMap[name] = val

        if(printMap): print(valMap)

        # set variables
        for obj in self.flexures + self.stages:
            obj.setPrms(valMap)

    def _resetVariables(self):
        
        for obj in self.flexures + self.stages:
            obj.resetPrms()

class CompMechSimp(object):
    def __init__(self, modelInfo, viewCenter):
        
        self.material = Material.create(["Isotropic", 0, "Dummy", 0, 0, 0, 0]) # dummy
        self.center = str2num(viewCenter)
        self.ready = False

        # flexures
        self.stages = []
        self.modes = []
        self.rods = []
        self.cables = []
        self.sensors = []

        self.pivots = []

        # spaces
        self.modeFreedomSpaces = []
        self.modeConstraintSpaces = []
        self.modeAxes = []
        self.freedomUnion = None
        self.consIntersection = None
        self.rodSpace = None
        self.cableSpaces = []
        self.sensorSpaces = []

        # completion check
        self.rodComplete = False
        self.cableComplete = False

        # systemBases
        self.modesBases = []
        self.rodsBases = []
        self.cablesBases = []
        self.sensorsBases = []

        # messages
        self.msgModeOverall = ''
        self.msgMode = []
        self.msgRod = ''
        self.msgCableConfig = ''
        self.msgCable = []
        self.sensorResponse = ''
        self.msgSensor = []

        self.msgRodFinal = ''
        self.msgCableFinal = ''
        self.msgSensorFinal = ''

        self._model(modelInfo)

    def _model(self, modelInfo):

        modelInfo = str2num(modelInfo)
        
        if(len(modelInfo[0]) > 0):
            self._modelStages(modelInfo[0])
        
        # check input status
        modeCleaned = []
        for i in range(len(modelInfo[1])):
            mode = modelInfo[1][i]
            if(len(mode) > 0):
                modeCleaned += [mode]
        
        dofCount = 0
        for mode in modelInfo[1]:
            dofCount += len(mode)
        
        inputValid = dofCount != 0
        
        if(inputValid):
            self._modelModes(modeCleaned)
            self._modelRods(modelInfo[2])
            self._modelCables(modelInfo[3])
            self._modelSensors(modelInfo[4])
            self.ready = True
        
    def _modelStages(self, info):

        for stageDesign in info:
            stage = Stage.model(stageDesign, self.material)
            self.stages += [stage]
        
    def _modelModes(self, info):
        
        # get DOF vectors
        for mode in info:
            vectors = np.stack([lineToScrewVec(tm) for tm in mode])
            self.modes += [vectors]
        
        # covert into freedom spaces
        for i in range(len(self.modes)):
            mode = self.modes[i]
            freedomSpace = Subspace(spans=mode)
            constraintSpace = Subspace(spans=DOFSystem.freeToCon(freedomSpace.spanSpace))
            self.modeFreedomSpaces += [freedomSpace]
            self.modeConstraintSpaces += [constraintSpace]
            
            # find DOF axes
            sys = DOFSystem(freedomSpace.spanSpace)
            
            orthoDOFs = sys.orthoDOF()
            
            axes = np.stack([getAxis(vec, True) for vec in orthoDOFs])
            self.modeAxes += [axes]
            
        # find DOF union and DOC intersection
        allFreedomVecs = np.concatenate([freedomVecs for freedomVecs in self.modes])
        self.freedomUnion = Subspace(spans=allFreedomVecs)
        if(self.freedomUnion.spanRank == MAXDOFS):
            self.consIntersection = Subspace(normals=np.identity(MAXDOFS))
        else:
            self.consIntersection = Subspace(spans=DOFSystem.freeToCon(self.freedomUnion.spanSpace))
        
    def _modelRods(self, info):
        
        # create flexures
        self.rods = []
        for rodDesign in info:
            rod = Flexure.model(rodDesign, self.material)
            self.rods += [rod]
        
        # convert into constraint space
        if(len(self.rods) > 0):
            screwVecs = np.stack([rod.screwSpace() for rod in self.rods])
            conSpace = Subspace(spans=screwVecs)
            self.rodSpace = conSpace
        else:
            self.rodSpace = Subspace(np.identity(MAXDOFS))

    def _modelCables(self, info):
        
        # create flexures
        for group in info:
            temp = []
            for cableDesign in group:
                cable = Flexure.model(cableDesign, self.material)
                temp += [cable]
            self.cables += [temp]
        
        # create constarint space
        for group in self.cables:
            if(len(group) > 0):
                screwVecs = np.stack([cable.screwSpace() for cable in group])
                conSpace = Subspace(spans=screwVecs)
            else:
                conSpace = Subspace(np.identity(MAXDOFS))
            
            self.cableSpaces += [conSpace]
            
    def _modelSensors(self, info):
        
        # create flexures
        for group in info:
            temp = []
            for sensorDesign in group:
                sensor = Flexure.model(sensorDesign, self.material)
                temp += [sensor]
            self.sensors += [temp]
            
        # create constraint space
        for group in self.sensors:
            if(len(group) > 0): 
                screwVecs = np.stack([sensor.screwSpace() for sensor in group])
                conSpace = Subspace(spans=screwVecs)
            else:
                conSpace = Subspace(np.identity(MAXDOFS))
                
            self.sensorSpaces += [conSpace]
        
    def printInfo(self):

        # base info
        print(self.material)
        print(self.center)
        
        # elements info
        print(self.stages)
        print(self.modes)
        print(self.rods)
        print(self.cables)
        print(self.sensors)
        
        # spaces
        print(self.modeFreedomSpaces)
        print(self.modeConstraintSpaces)
        print(self.freedomUnion)
        print(self.consIntersection)
        print(self.rodSpace)
        print(self.cableSpaces)
        print(self.sensorSpaces)

    def outputMsg(self):

        if(not self.ready):
            return []

        msgs = []
        msgs += [tuple([self.msgModeOverall, tuple(self.msgMode)])]
        msgs += [tuple([self.msgRod])]
        msgs += [tuple([self.msgCableConfig, tuple(self.msgCable)])]
        msgs += [tuple([self.sensorResponse, tuple(self.msgSensor)])]
        msgs += [tuple([self.msgRodFinal, self.msgCableFinal, self.msgSensorFinal])]

        return tuple(msgs)

    def outputBases(self):
        
        output = []
        output += [tuple(self.modesBases)]
        output += [self.rodsBases]
        output += [tuple(self.cablesBases)]
        output += [tuple(self.sensorsBases)]

        return tuple(output)

    def computeDesign(self):

        if(not self.ready):
            return []

        feedback = []
        # for step 1 freedom space visualization
        feedback += [self._computeModeFreedomSpaces()]
        # for step 2 rod constraint space visualization
        feedback += [self._computeConsIntersection()]
        # for step 3 cable modeling - constraint space visualization
        feedback += [self._computeCableGroups()]
        # for step 4 sensor modeling
        feedback += [self._computeSensorGroups()]

        self._completionCheck()
        
        return tuple(feedback)

    def _computeModeFreedomSpaces(self):
        # overall viability check
        allFreedomSys = DOFSystem(self.freedomUnion.spanSpace)
        self.msgModeOverall = Msngr.viabilityCheck(allFreedomSys)

        spaces = []
        # output freedom spaces for each
        for mode in self.modeFreedomSpaces:
            sys = DOFSystem(mode.spanSpace)
            
            if(self.center != None):
                sys.recenterSpaces(self.center)
            
            
            motions = tuple(sys.output())
            self.modesBases += [sys.outputBases()]

            spaces += [motions]

            msg = Msngr.describeFreedomSpace(sys, motions)
            self.msgMode += [msg]

        return tuple(spaces)

    def _computeConsIntersection(self): 
        
        msg = [Msngr.showNeeded]

        # check overconstraints
        msgOC, ocRods = self._rodOverConstraint()
        noOC = len(ocRods) == 0

        # check completion
        msgComp = self._rodCompletion()

        # show completion if not over-constrained
        if(noOC):
            msg += [msgComp]
        else:
            msg += msgOC
        
        # get constraint space
        sys = DOFSystem(self.consIntersection.spanSpace, True)

        if(self.center != None):
            sys.recenterSpaces(self.center)
        constraints = tuple(sys.output())
        
        self.rodsBases = sys.outputBases()

        # describe space
        msg += Msngr.describeRodSpace(constraints)
        msg = '\n'.join(msg)
        
        self.msgRod = msg

        return tuple([ocRods, constraints])

    def _rodCompletion(self):

        tarSpace = self.consIntersection.spanSpace
        curSpace = self.rodSpace.spanSpace
        
        tarSys = DOFSystem(tarSpace, True)
        curSys = DOFSystem(curSpace, True)

        tarAxisDeg, tarSpanDeg, tarRods = tarSys.rotDeg, tarSys.transDeg, tarSys.deg
        curAxisDeg, curSpanDeg, curRods = curSys.rotDeg, curSys.transDeg, curSys.deg

        msg = Msngr.completion(tarAxisDeg, tarSpanDeg, tarRods, curAxisDeg, curSpanDeg, curRods)

        completion = tarAxisDeg == curAxisDeg and\
                     tarSpanDeg == curSpanDeg and\
                     tarRods == curRods
        
        status = Msngr.rodSpaceComplete if completion else Msngr.rodSpaceIncomplete
        msg = '\n'.join([status, msg])

        self.rodComplete = completion
        
        return msg

    def _rodOverConstraint(self):

        targSpace = self.consIntersection.spanSpace
        msg = []

        # check calidity of each rod
        isValid = True
        invalidRods = []
        for rod in self.rods:
            isAllowed = rod.isAllowed(targSpace)
            if(not isAllowed):
                invalidRods += [rod]
                isValid = False
        
        if(isValid):
            msg += [Msngr.noOCRod]
        else:
            msg += [Msngr.hasOCRod]
            msg += Msngr.ocIssue(invalidRods, targSpace, "Rod")
        
        ocRods = tuple([rod.id for rod in invalidRods])
        
        return msg, ocRods

    def _computeCableGroups(self):

        # get cable groups
        subspaces = self._computeConstraintSubspaces() # narrow
        
        # find group's twin pivot axes
        pivots = self._cablePivotAxes(subspaces)

        # generate reconfiguration plan
        self.msgCableConfig = self._genReconfigPlan(subspaces)
        
        # find cable placement spaces
        groupSpaces = self._expandedCableSpace(subspaces) # extended
        
        # check over-constraints
        msgOC, inValidCables = self._cableOverConstraint(groupSpaces)
        noOC = [len(ocCables) == 0 for ocCables in inValidCables]
        
        # check completion
        msgComp, highlights = self._cableCompletion(subspaces, groupSpaces)

        # create and decribe constraint spaces
        conSpaces, msgCon = self._cableSpaceOutput(groupSpaces)

        # compose messages
        msgs = []
        for i in range(len(msgOC)):
            if(noOC[i]):
                groupMsgs = [msgOC[i]] + msgComp[i] + msgCon[i]
            else:
                groupMsgs = [msgOC[i]]
            msg = '\n'.join(groupMsgs)
            msgs +=[msg]
        
        self.msgCable = msgs

        return tuple([inValidCables, highlights, conSpaces, pivots])

    def _computeConstraintSubspaces(self):

        # find all combinations of modes
        modeCount = len(self.modes)
        ids = range(modeCount)
        combs = [] # sorted by level
        for i in range(modeCount):
            lvlCombs = combinations(ids, i + 1)
            combs += [comb for comb in lvlCombs]
        
        # intersect all combinations to find all subspaces in the Venn diagram
        subspaces = {}
        for comb in combs:
            freedomSpaces = [self.modes[id] for id in comb]
            freedomConcat = np.concatenate(freedomSpaces)
            freedomAll = Subspace(spans=freedomConcat)
            conSpace = DOFSystem.freeToCon(freedomAll.spanSpace)\
            
            if(conSpace.size == 0):
                space = Subspace(normals=np.identity(MAXDOFS))
            else:
                space = Subspace(spans=conSpace)
            
            subspaces[comb] = space
        
        # simplify spaces by removing subset spaces from a space
        # if a space ends up empty, then it's removed to simplify the design
        emptySubspaces = []
        for comb in reversed(subspaces): # starting from the highest level
            
            # find all subsets of a constraint space
            subset = [x for x in combs[combs.index(comb) + 1:]\
                         if len(x) > len(comb) and\
                            set(comb).issubset(set(x))]
            
            if(len(subset) == 0): continue # no subset, skip
            
            # gather all subset constraint screw vectors and simplify
            subsetSpaces = [subspaces[subsetKey].spanSpace for subsetKey in subset]
            subsetSpaces = [space for space in subsetSpaces if space.size != 0]
            
            if(len(subsetSpaces) != 0):
                subsetConcat = np.concatenate(subsetSpaces)
                subsetSpace = Subspace(spans=subsetConcat)
                subsEmpty = False
            else:
                subsetSpace = Subspace(normals=np.identity(MAXDOFS))
                subsEmpty = True
            
            # compute the difference
            tarSpace = subspaces[comb]
            if(subsEmpty):
                diff = tarSpace
            else:
                diff = Subspace.difference(tarSpace, subsetSpace)
            diff.subsetSpaces = [subspaces[sub] for sub in subset]
            subspaces[comb] = diff

            # flag the subspace if it ends up empty for later removal
            if(diff.spanRank == 0):
                emptySubspaces += [comb]
        
        # remove empty subspaces
        for comb in emptySubspaces:
            subspaces.pop(comb, None)
        
        for comb in subspaces:
            subspace = subspaces[comb]
            subspace.subsetSpaces = [space for space in subspace.subsetSpaces if space.spanRank > 0]

        # pop all-intersection(i.e., the rods)
        subspaces.pop(tuple(range(modeCount)), None)

        return subspaces
    
    def _genReconfigPlan(self, subspaces):

        modeCount = len(self.modes)
        groupCount = len(subspaces)

        rowTitle = ' ' * 5 + ''.join([" Group%d" % (i + 1) for i in range(groupCount)])
        msg = [rowTitle]
        for modeId in range(modeCount):
            row = "Mode%d" % (modeId + 1)
            for subspaceIds in subspaces:
                if(modeId in subspaceIds):
                    row += "  Tight"
                else:
                    row += "  Loose"
            msg += [row]
        
        return '\n'.join(msg)

    def _expandedCableSpace(self, subspaces):

        groupSpaces = {}
        for group in subspaces:
            curSpace = subspaces[group]
            subsetAll = [space.spanSpace for space in curSpace.subsetSpaces]

            if(len(subsetAll) == 0):
                groupSpaces[group] = curSpace
                continue

            subsetAll = np.concatenate(subsetAll)
            subsetSpace = Subspace(spans=subsetAll)
            unionSpace = Subspace(spans=np.concatenate([curSpace.spanSpace, subsetSpace.spanSpace]))
            unionSys = DOFSystem(unionSpace.spanSpace, True)
            freeRots = unionSys.findFreeRotations()
            groupSpace = Subspace(spans=np.concatenate([curSpace.spanSpace, freeRots]))
            
            groupSpaces[group] = groupSpace

        return groupSpaces

    def _cableOverConstraint(self, groupSpaces):
        
        msgs = []
        inValidCables = []

        # bypass
        if(len(self.cables) > len(groupSpaces.keys())):
            return msgs, inValidCables
        
        # check overConstraint
        for i in range(len(self.cables)):
            isValid = True
            temp = []
            key = list(groupSpaces.keys())[i]

            # create twins

            for j in range(len(self.cables[i])):
                cable = self.cables[i][j]
                expandedSpace = groupSpaces[key].spanSpace
                
                isAllowed = cable.isAllowed(expandedSpace)
                if(not isAllowed):
                    temp += [cable]
                    isValid = False
            
            inValidCables += [tuple([cable.id for cable in temp])]

            msg = []            
            if(isValid):
                msg += [Msngr.noOCCable]
            else:
                msg += [Msngr.hasOCCable]
                msg += Msngr.ocIssue(temp, groupSpaces[key].spanSpace, "Cable")

            msgs += ['\n'.join(msg)]


        
        return msgs, inValidCables

    def _cableCompletion(self, subspaces, groupSpaces):

        msgs = []
        highlights = []
        done = []
        subspacesKeys = list(subspaces.keys())
        
        for i in range(len(subspacesKeys)): # *needed* only
            # target space
            tarSpace = subspaces[subspacesKeys[i]]
            
            
            # gather all subset constraint spaces
            allVecs = [space.spanSpace for space in tarSpace.subsetSpaces]
            # add cables if present
            if(i < len(self.cables)):
                # find the cable space with twins
                cableSpace = self._cableTwinCompleteSpace(i)
                #cableSpace = self.cableSpaces[i]
                
                if(cableSpace != None and cableSpace.spanRank != 0):
                    allVecs += [self.cableSpaces[i].spanSpace]
            
            # calculate the difference
            if(len(allVecs) != 0):
                allVecs = np.concatenate(allVecs)
                curSpace = Subspace(spans=allVecs)
                diff = Subspace.difference(tarSpace, curSpace)
            else:
                diff = tarSpace
            
            isComplete = diff.spanRank == 0
            if(isComplete):
                # show completion message
                msg = [Msngr.cableSpaceComplete]
                highlight = []
            else:
                # generate highlight and decribe missing part
                diffSys = DOFSystem(diff.spanSpace)
                msg = [Msngr.neededCable(diffSys.rotDeg, diffSys.transDeg)]
                highlight, msgPerp = self._getHighlights(diff)
                if(len(highlight) != 0): msg += [msgPerp]
                msg += [Msngr.cableSpaceInComplete]

            msgs += [msg]
            done += [isComplete]
            highlights += [highlight]
        
        self.cableComplete = False not in done
        return msgs, highlights

    def _cableTwinCompleteSpace(self, i):

        cables = self.cables[i]
        axes = self.pivots[i]

        center = np.asarray(self.center) if self.center != None else np.zeros(EUCSPACEDIM)

        stack = []
        for cable in cables:
            start, end = cable.start, cable.end
            startAll, endAll = [start], [end]
            for axis in axes:
                
                rotSV = screwVec(axis, center) * np.pi
                
                startNew, endNew = [], []
                for i in range(len(startAll)):
                    twinStart = screwTransform(startAll[i], rotSV)
                    twinEnd = screwTransform(endAll[i], rotSV)
                    startNew += [twinStart]
                    endNew += [twinEnd]
                
                startAll += startNew
                endAll += endNew
            
            startAll = np.asarray(startAll)
            endAll = np.asarray(endAll)
            vecAll = endAll - startAll
            stack += [screwVec(vecAll[i], startAll[i]) for i in range(len(startAll))]
        
        if(len(stack) == 0):
            return None
        else:
            fullSpace = Subspace(spans=np.stack(stack))
            
            return fullSpace

    def _getHighlights(self, diff):
        
        diffSys = DOFSystem(diff.spanSpace, True)

        if(diffSys.rotDeg == 0 or diffSys.transDeg > 0): # don't need unique directions, skip
            return (), ""

        # find all needed new directions        
        neededDir = []
        for rotDOF in diffSys.permRots:
            neededDir += [list(rotDOF[1])]
        neededDir = np.stack(neededDir)
        
        neededDir = Subspace(spans=neededDir).spanSpace
        neededDir = neededDir / np.linalg.norm(neededDir, 1)
        neededDir = tuple(map(tuple, neededDir))
        
        return neededDir, Msngr.cableNotPerp

    def _cableSpaceOutput(self, groupSpace):

        spaces = []
        msgs = []
        for group in groupSpace:
            
            sys = DOFSystem(groupSpace[group].spanSpace, True)
            
            if(self.center != None):
                sys.recenterSpaces(self.center)
        
            constraints = tuple(sys.output())
            
            self.cablesBases += [sys.outputBases()]

            msg = Msngr.describeCableSpace(constraints)

            spaces += [constraints]
            msgs += [msg]


        return spaces, msgs

    def _cablePivotAxes(self, subspaces):
        
        allModes = range(len(self.modes))
        
        groups = subspaces.keys()

        pivots = []
        for group in groups:
            temp = []
            for id in allModes:
                if(id in group): continue

                axes = self.modeAxes[id]
                temp += [axes]
            
            temp = np.concatenate(temp)
            temp = Subspace(spans=temp).spanSpace
            
            self.pivots += [temp]
            
            axes = [tuple(vec) for vec in temp]
            pivots += [tuple(axes)]

        return tuple(pivots)

    def _computeSensorGroups(self):

        # find cable placement spaces
        groupSpaces = self._expandedSensorSpace() # extended

        # find mode's twin pivot axes
        pivots = self._sensorPivotAxes()
        
        # check responsiveness
        responsiveness = self._responsiveness()
        self.sensorResponse = Msngr.responseMsg(responsiveness)
        
        # create and decribe constraint spaces
        conSpaces, msgCon = self._sensorSpaceOutput(groupSpaces)
        
        # sensor icon info
        icon = self._sensorIcon(responsiveness)

        # compose messages
        for i in range(len(groupSpaces)):
            header = Msngr.sensorShow(i)
            description = msgCon[i]
            msg = '\n'.join([header] + description)
            self.msgSensor += [msg]

        return tuple([pivots, conSpaces, icon])

    def _expandedSensorSpace(self):

        spaces = []
        for i in range(len(self.modes)):
            mode = self.modeFreedomSpaces[i]
            # get from mode's invalid constraint space
            conSpace = Subspace(spans=DOFSystem.freeToCon(mode.spanSpace))
            conSys = DOFSystem(conSpace.spanSpace, True)
            invalidConSpace = conSpace.normalSpace
            
            # supplement directional vectors from mode's constraint space
            center = self.center if self.center != None else np.zeros(EUCSPACEDIM)
            freeRots = conSys.basesRots(center)
            
            allVecs = np.concatenate([invalidConSpace, freeRots])
            allVecs = Subspace(spans=allVecs).spanSpace
            spaces += [allVecs]
            
        return spaces

    def _sensorPivotAxes(self):

        allAxes = []
        for mode in self.modeAxes:
            temp = []
            for axis in mode:
                temp += [tuple(axis)]
            allAxes += [tuple(temp)]

        return tuple(allAxes)
        
    def _responsiveness(self):
        
        modeFreedoms = []

        responsiveness = []
        # for each kinematic mode
        for i in range(len(self.modes)):
            mode = self.modeFreedomSpaces[i]
            modeSys = DOFSystem(mode.spanSpace)
            freedoms = modeSys.indiDOFs()
            modeFreedoms += [freedoms]

            resMap = {}
            for dof in freedoms:
                
                freedom = freedoms[dof]
                constraint = DOFSystem.freeToCon(freedom)
                conRank  = len(constraint)

                temp = []
                # for each cable group
                for j in range(len(self.sensorSpaces)):
                    sensorSpace = self.sensorSpaces[j].spanSpace
                    responsive = self._checkResponsive(sensorSpace, constraint)
                    temp += [responsive]
                
                resMap[dof] = temp

            responsiveness += [resMap]
        
        return responsiveness
    
    def _checkResponsive(self, sensorSpace, constraint):

        conRank  = len(constraint)
        if(sensorSpace.size == 0): # no sensor
                return False
        
        # if the sensor space is not allowed by the freedom space's
        # corresponding constraint space, then it is responsive

        # checking subsets by combined subspace rank increase (or not)
        combined = np.concatenate([sensorSpace, constraint])
        combinedRank = Subspace(spans=combined).spanRank
        responsive = combinedRank > conRank
        
        return responsive

    def _sensorIcon(self, responsiveness):

        modeCount = len(responsiveness)

        # generate container for each mode
        iconStacks = []
        for i in range(modeCount):
            iconStacks += [["base.png"]]

        # put dofs into individual modes
        for i in range(modeCount):
            modeDofs = responsiveness[i]
            for dof in modeDofs:
                iconStacks[i] += [dof + ".png"]
                if(modeDofs[dof][i]):
                    iconStacks[i] += [dof + '_on.png']

        iconStacks = [tuple(iconStacks[i]) for i in range(len(iconStacks))]
        
        return tuple(iconStacks)

    def _sensorSpaceOutput(self, groupSpaces):
        
        spaces = []
        msgs = []
        for space in groupSpaces:
            sys = DOFSystem(space, True)

            if(self.center != None):
                sys.recenterSpaces(self.center)
        
            conSpaces = tuple(sys.output())
            
            self.sensorsBases += [sys.outputBases()]

            msg = Msngr.describeSensorSpace(conSpaces)
            
            spaces += [conSpaces]
            msgs += [msg]
        # for each kinematic mode and sensor group
        # render space and generate description
        
        return tuple(spaces), msgs

    def _completionCheck(self):

        # rods
        if(self.rodComplete):
            msg = [Msngr.rodComplete] + Msngr.rodRecommendDia(self.rods)
            self.msgRodFinal = '\n'.join(msg)
        else:
            self.msgRodFinal = Msngr.rodIncomplete
        
        # cables
        if(self.cableComplete):
            msg = [Msngr.cableComplete, self.msgCableConfig]
            self.msgCableFinal = '\n'.join(msg)
        else:
            self.msgCableFinal = Msngr.cableIncomplete

        # sensors
        self.msgSensorFinal = self.sensorResponse

class Msngr(object):
    showNeeded = MSG_SHOW_NEEDED
    noOCRod = MSG_NO_OC_ROD
    hasOCRod = MSG_HAS_OC_ROD
    noOCCable = MSG_NO_OC_CABLE
    hasOCCable = MSG_HAS_OC_CABLE
    rodComplete = MSG_ROD_COMPLETE
    rodIncomplete = MSG_ROD_INCOMPLETE
    rodSpaceComplete = MSG_ROD_SPACE_COMPLETE
    rodSpaceIncomplete = MSG_ROD_SPACE_INCOMPLETE

    cableComplete = MSG_CABLE_COMPLETE
    cableIncomplete = MSG_CABLE_INCOMPLETE
    cableSpaceComplete = MSG_CABLE_SPACE_COMPLETE
    cableSpaceInComplete = MSG_CABLE_SPACE_INCOMPLETE
    cableNotPerp = MSG_CABLE_NOT_PERP

    sensorSuggest = MSG_SENSOR_SHOW
    sensorNote = MSG_SENSOR_NOTE

    def __init__(self):
        pass

    @staticmethod
    def describeSubspace(motion, name=None, isFreedom=True):
        
        moType = motion[0]
        axisDeg = 0
        for axis in motion[1]:
            if(len(axis) == 3):
                axisDeg += 1
        spaceDeg = 0
        for axis in motion[3]:
            if(len(axis) == 3):
                spaceDeg += 1
        
        # space type
        mode = None
        if(moType == 'r' and isFreedom):
            spaceType = "rotation"
            mode = 'r'
        elif(moType == 't' and isFreedom):
            spaceType = "translation"
            axisDeg = spaceDeg
            mode = 't'
        elif(moType == 'r' and not isFreedom):
            spaceType = "wrench/wire"
            mode = 'w'
        
        # name override
        if(name != None):
            spaceType = name
        
        # composition
        axisNote = ''
        if(axisDeg == 1):
            direction = "that aligns with the displayed direction"
        elif(axisDeg == 2):
            direction = "that is parallel to the plane(s)"
        elif(axisDeg == 3):
            direction = "in 3D space"
            axisNote = "Note: axes must not all lie on the same plane."
        
        spaceNote = ''
        # location):
        if(spaceDeg == 0):
            position = "passes through the center point"
        elif(spaceDeg == 1):
            position = "passes through this line at some point"
        elif(spaceDeg == 2):
            if(axisDeg <= 2):
                position = "lies on this plane"
            else:
                position = "intersects with this plane at some point"
            if(axisDeg == 2 and mode != 't'):
                spaceNote = "Note: the axes must not intersect at the same point."
        elif(spaceDeg == 3):
            if(axisDeg == 1):
                position = "lies anywhere in space"
                if(mode != 't'):
                    spaceNote = "Note: axes must not all lie on a same plane."
            elif(axisDeg == 2):
                position = "lies on any parallel plane"
                if(mode != 't'):
                    spaceNote = "Note: axes must not all lie on the same plane."
            elif(axisDeg == 3):
                position = "lies anywhere in space"
                if(mode != 't'):
                    spaceNote = "Note: axes must not all lie on a same plane nor pass through the same point."
        
        if(moType == 'r'):
            msg = "Any axis %s and %s is an allowed %s axis. " % (direction, position, spaceType)
        elif(moType == 't'):
            msg = "Any axis %s is an allowed %s axis. " % (direction, spaceType)
        
        
        msg += axisNote + ' ' + spaceNote
        
        return msg

    @staticmethod
    def describeFreedomSpace(sys,motions):
        
        msgs = []
        # report number of freedoms
        msgs += ["This mode has %d rotational and %d translational DOF." % (sys.rotDeg, sys.transDeg)]
        
        # check viability
        msgs += [Msngr.viabilityCheck(sys, "Current mode")]
        
        
        # describe the direction and location of subspaces
        for i in range(len(motions)):
            motion = motions[i]
            msgs += ["Motion subspace %d:" % (i + 1)]
            msgs += [Msngr.describeSubspace(motion)]

        return '\n'.join(msgs)

    @staticmethod
    def describeRodSpace(motions):
        
        # describe direction and position
        msg = []

        for i in range(len(motions)):
            con = motions[i]
            msg += ["Constraint subsapce %d:" % (i + 1)]
            msg += [Msngr.describeSubspace(con, "rod", False)]
        
        return msg

    @staticmethod
    def describeCableSpace(motions):
        
        # describe direction and position
        
        msg = []

        for i in range(len(motions)):
            con = motions[i]
            msg += ["Constraint subsapce %d:" % (i + 1)]
            msg += [Msngr.describeSubspace(con, "cable", False)]
        
        return msg

    @staticmethod
    def describeSensorSpace(motions):
        
        # describe direction and position
        
        msg = []

        for i in range(len(motions)):
            con = motions[i]
            msg += ["Constraint subsapce %d:" % (i + 1)]
            msg += [Msngr.describeSubspace(con, "sensor", False)]
        
        return msg
    
    @staticmethod
    def viabilityCheck(sys, name="Design"):
        
        if(sys.isConstraint):
            isViable = sys.rotDeg > 0 # has directional component
        else:
            # is not unconstrained & has directional component
            isViable = (sys.rotDeg + sys.transDeg < 6) and (sys.transDeg < 3)
        
        viability = "viable" if isViable else "unviable"
        msg = "%s is %s." % (name, viability)

        return msg

    @staticmethod
    def ocIssue(elems, space, name):

        msgs = []
        for elem in elems:
            validDir, validPos = elem.dirPosAllowed(space)
            if(not validDir and validPos):
                violation = "direction"
            elif(validDir and validPos):
                violation = "position"
            else:
                violation = "direction and position"

            msg = "The highlighted %s has an invalid %s." % (name + " %d" % elem.id, violation)
            msgs += [msg]
        
        return msgs

    @staticmethod
    def completion(*args):

        return MSG_CHECKLIST(*args)

    @staticmethod
    def responseMsg(responsiveness):

        modeCount = len(responsiveness)
        
        # generate matrix
        resMatrix = []
        for i in range(modeCount):
            modeDofs = responsiveness[i]
            res = [False] * modeCount
            
            for dof in modeDofs:
                resDof = modeDofs[dof]
                for j in range(modeCount):
                    if(resDof[j]):
                        res[j] = True
            
            resMatrix += [res]
        
        # generate message
        rowTitle = ' ' * 5 + ''.join([" Group%d" % (i + 1) for i in range(modeCount)])
        msg = [rowTitle]
        for modeId in range(modeCount):
            row = "Mode%d" % (modeId + 1)
            for groupId in range(modeCount):
                responsive = resMatrix[modeId][groupId]
                if(responsive):
                    row += "   V   "
                else:
                    row += "   X   "
            msg += [row]
        
        msg += [Msngr.sensorNote]
        
        return '\n'.join(msg)

    @staticmethod
    def sensorShow(id):

        return Msngr.sensorSuggest + " Mode %d." % id

    @staticmethod
    def rodRecommendDia(rods):
        
        if(len(rods) == 0):
            avgLen = 0
        else:
            lengths = []
            for rod in rods:
                lengths += [rod.length]
            
            avgLen = sum(lengths) / len(lengths)

        minLen = avgLen / ROD_MIN_ASPECT_RATIO
        maxLen = avgLen / ROD_MAX_ASPECT_RATIO

        msg = ["Average flexural rod length: %.2f mm." % avgLen]
        msg += ["Recommended rod diameter: %.2f - %.2f mm" % (minLen, maxLen)]
        msg += ["You can extend or trim the curves without changing the device's kinematics."]
        
        return msg
    
    @staticmethod
    def neededCable(rotDeg, transDeg):
        
        if(rotDeg > 0 and transDeg == 0):
            needed = "%d more unique directions" % rotDeg
        elif(rotDeg == 0 and transDeg > 0):
            needed = "%d more unique positions" % transDeg
        elif(rotDeg > 0 and transDeg > 0):
            needed = "%d and %d more unique directions and positions, respectively" % (rotDeg, transDeg)
        # the other one is an impossible case: nothing needed
        
        msg = "Need %s." % needed

        return msg

if __name__ == "__main__":

    m = CompMechSimp(inputModel, viewCenter)
    info = m.computeDesign()
    messages = m.outputMsg()
    bases = m.outputBases()