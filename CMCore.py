"""Provides a scripting component.
    Inputs:
        x: The x script variable
        y: The y script variable
    Output:
        a: The a output variable"""

__author__ = "MorphingMatterLab"
__version__ = "2021.03.16"

SCALE = 10
EPSILON = 1e-6
DEBUG_MODE = False
MSG_UNKNOWN_CON_TYPE = "Unsupported object found during constraint modeling"
MSG_UNKNOWN_STAGE_TYPE = "Unknown stage type"
MSG_UNIMPLEMENTED = "Funciton not implemented"

DEFAULT_WIRE_RAD = 1

import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino
import math
from time import time
from decimal import *

SimpMode = True
getcontext().prec = 6
################################################################################
# base classes
################################################################################

class InputObject(object):
    def __init__(self, objType, geo=None, adj=None):
        self.geo = geo
        self.adj = adj
        self.type = objType
        self.deleteFlag = False
    
    def __eq__(self, other):
        
        if(isinstance(other, InputObject)):
            return self.type == other.type and self.geo == other.geo
        if(isinstance(other, tuple)):
            return self.adj == other
        else:
            return self.geo == other
            
    def checkValidity(self):
        
        if(self.deleteFlag): return False, 0
        
        if(self.adj != None): # if not stage
            if(isinstance(self.adj, tuple)): # flexure elements
                for obj in self.adj:
                    isValid = obj.checkValidity()
                    if(not isValid): return False, 2
            else: # goals
                isValid = self.adj.checkValidity()
                if(not isValid): return False, 2
        
        # stage check self or goals, flexures check self
        
        if(self.geo == None):
            isValid = True
        else:
            isValid = InputObject._checkObject(self.geo)
        
        if(isValid): return True, None
        else: return False, 1
    
    def __repr__(self):

        if(self.adj == None):
            msg = "InputObject(%s, %s)" % (self.type, self.geo)
        else:
            msg = "InputObject(%s, %s, %s)" % (self.type, self.geo, self.adj)
        
        return msg

    @staticmethod
    def _checkObject(obj):
        
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        
        # check if object exists
        o = rs.coercegeometry(obj)
        if(o == None):
            isValid = False
        else:
            isValid = True
        
        sc.doc = ghdoc
        
        return isValid

class Constraint(object):
    name = "Constraint"
    degree = 0
    isElemental = True
    types = []
    initialized = False
    
    def __init__(self, geo, id, adjs=None):
        
        self.geo = geo
        self.id = id
        self.adjs = adjs
        self.mat = 0 # hardcoded
        
        if(not Constraint.initialized):
            Constraint._initConstraintTypes()
    
    @staticmethod
    def _initConstraintTypes():
         
        for subclass in Constraint.__subclasses__():
            if(not subclass.isElemental): continue
            Constraint.types += [subclass]
        
        Constraint.initialized = True
    
    def __repr__(self):
        
        msg = "<Constraint #%d" % self.id
        if(DEBUG_MODE):
            msg += " adj:%s" % str(self.adjs)
        
        msg += '>'
        
        return msg
    
    def getGeometries(self):
        
        geometries = self.geo
        
        if(not isinstance(geometries, list)):
            geometries = [geometries]
        
        return geometries
    
    def findAdjStages(self, stages):
        
        adjId = []
        for id in stages:
            stage = stages[id]
            for solid in stage.solids:
                if(self._isAdjacent(solid.geo)):
                    adjId += [stage.id]
                    break
        
        adjId = tuple(sorted(adjId))
        
        return adjId
    
    def _isAdjacent(self, solid):
        # place holder function
        # inherit in child classes
        
        assert False, MSG_UNIMPLEMENTED
        
        return 42
    
    def freedom(self):
        # place holder function
        # inherit in child classes
        
        assert False, MSG_UNIMPLEMENTED
        
        return 42 # Matrix of size (N, 6)
    
    def outputBasicInfo(self):
        
        info = []
        info += [type(self).name]
        info += [self.id]
        info += [tuple([s.id for s in self.adjs])]
        info += [self.mat]
        
        return info
    
    def outputInfo(self):
        # place holder function
        # inherit in child classes
        
        assert False, MSG_UNIMPLEMENTED
        
        return 42 # tuple
    
    def outputTags(self):
        
        center = self._tagPt()
        
        id = self.id
        
        return center, id
    
    def outputVis(self):
        
        return [self.geo]
    
    def _tagPt(self):
        # place holder 
        # inherit in child classes
        
        assert False, MSG_UNIMPLEMENTED
        
        return rs.CreatePoint(0, 0, 0) # returns a Point3d
    
    @staticmethod
    def model(*args):
        
        converted = None
        typeName = args[0]
        subclasses = exhaustiveSubclasses(Constraint)

        for conType in subclasses:
            if(typeName == conType.name):
               converted = conType(*args[1:])
               break
        
        return converted
    
    @staticmethod
    def _isType(geometry):
        # place holder 
        # inherit in child classes
        
        assert False, MSG_UNIMPLEMENTED
        
        return False # returns a bool value

class Joint(object):
    name = "Joint"
    
    def __init__(self, cons, id, adjs):
        self.cons = cons
        self.id = id
        self.adjStages = adjs
    
    def __repr__(self):
        
        msg = "<Joint #%d" % self.id
        if(DEBUG_MODE):
            msg += " #C:%d adj:%s" % (len(self.cons),
                                      str([s.id for s in self.adjStages]))
        
        msg += '>'
        
        return msg
        
    def outputInfo(self):
        
        info = []
        info += [Joint.name]
        info += [self.id]
        info += [tuple([s.id for s in self.adjStages])]
        info += [tuple([c.id for c in self.cons])]
        
        return tuple(info)
    
    def outputTags(self):
        
        center = self._getCenter()
        id = self.id
        
        return center, id
    
    def _getCenter(self):
        
        p0 = self.adjStages[0].center
        p1 = self.adjStages[1].center
        center = rs.PointDivide(rs.PointAdd(p0, p1), 2)
        
        return center

class Stage(object):
    name = "RigidSolid"
    
    def __init__(self, geo, id, isFixed=False):
        
        isDummy = geo == None

        if(not isinstance(geo, list) and not isinstance(geo, tuple)):
            geo = [geo]
        
        self.geo = geo
        self.geoCount = len(geo)
        self.id = id
        self.mat = 0
        
        self.isFixed = isFixed
        
        self.center = None
        self.volume = None
        
        if(isDummy):
            self.center = (0, 0, 0)
            self.volume = 0
        else:
            self._initPt()
    
    def _initPt(self):
        # weighted average method of computing centroid
        
        geometries = self.getGeometries()
        
        pts, vols = [], []
        for geometry in geometries:
            if(rs.IsBrep(geometry)): # brep method
                volume = rs.SurfaceVolume(geometry)[0]
                centroid = rs.SurfaceVolumeCentroid(geometry)[0]
            elif(rs.IsMesh(geometry)): # mesh method
                volume = rs.MeshVolume(geometry)[1]
                centroid = rs.MeshVolumeCentroid(geometry)
            else:
                assert False, MSG_UNIMPLEMENTED
            
            pts += [centroid]
            vols += [volume]
            
        pt = rs.CreatePoint(0, 0, 0)
        volSum = 0
        for i in range(len(geometries)):
            weighted = rs.PointScale(pts[i], vols[i])
            pt = rs.PointAdd(pt, weighted)
            volSum += vols[i]
        
        
        self.center = rs.PointDivide(pt, volSum)
        self.volume = volSum
        
    def __repr__(self):
        
        msg = "<Stage #%d" % self.id
        if(DEBUG_MODE):
            msg += " #solids:%d" % self.geoCount
        
        msg += '>'
        
        return msg
    
    def getGeometries(self):
        
        return self.geo
    
    def intersect(self, line):
        
        stage = self.geo[0]
        if(rs.IsMesh(stage)):
            closestPt = rs.CurveMeshIntersection(line, stage)
            
        else:
            intersection = rs.CurveBrepIntersect(line, stage)
            
            if(intersection == None):
                return None
                
            midPt = rs.CurveMidPoint(line)
            closest, minDist = None, None
            
            for intObj in intersection[1]:
                if(not rs.IsPoint(intObj)):
                    continue
                dist = rs.Distance(intObj, midPt)
                if(closest == None or dist < minDist):
                    closest, minDist = intObj, dist
            
            closestPt = [coord for coord in rs.coerce3dpoint(closest)]
            
            rs.DeleteObjects(intersection[1])
            
            return closestPt
    
    def minDistTo(self, pt):
        
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        
        shape = self.geo[0]
        if(rs.IsMesh(shape)):
            result = rs.MeshClosestPoint(shape, pt)
            cp = result[0]
        else:
            result = rs.BrepClosestPoint(shape, pt)
            cp = result[0]
        
        dist = rs.Distance(pt, cp)
        
        sc.doc = ghdoc
        
        return dist
    
    def outputTags(self):
        
        pt = self.center
        id = self.id
        
        return pt, id
    
    def outputInfo(self):
        
        info = []
        info += [Stage.name]
        info += [self.id]
        info += [self.mat]
        info += [p2t(self.center)]
        info += [str(dFix(self.volume))]
        info += [1 if self.isFixed else 0]
        
        return tuple(info)
    
    def outputVis(self):
        
        return self.geo
    
    @staticmethod
    def model(typeName, geo, id):
        if(typeName == Stage.name):
            return Stage(geo, id)
        else:
            return None

class Material(object):
    def __init__(self, id, name):
        self.id = id
        self.name = name
    
    def __eq__(self, other):
        
        return self.name == other.name

class IsoMat(Material):
    name = "Isotropic"
    
    def __init__(self, name, id, G, E, D, maxStrain):
        
        super(IsoMat, self).__init__(id, name)
        
        self.G = G # shear modulus
        self.E = E # elastic modulus
        self.D = D # density
        self.maxStrain = maxStrain # max strain

    def outputInfo(self):
        
        info = []
        info += [IsoMat.name]
        info += [self.id]
        info += [self.name]
        info += [str(dFix(self.E))]
        info += [str(dFix(self.G))]
        info += [str(dFix(self.D))]
        info += [str(dFix(self.maxStrain))]
        
        return tuple(info)

class Goal(object):
    def __init__(self, axis, id, stage):
        self.type = None
        self.id = id
        self.start = rs.CurveStartPoint(axis)
        self.end = rs.CurveEndPoint(axis)
        self.stage = stage
    
    def axis(self):
        
        return rs.AddLine(self.start, self.end)
    
    @staticmethod
    def model(*args):
        
        converted = None
        typeName = args[0]
        
        for subType in Goal.__subclasses__():
            if(typeName == subType.name):
               converted = subType(*args[1:])
               break
        
        return converted

    def outputInfo(self):
        
        info = []
        info += [self.type]
        info += [p2t(self.start)]
        info += [p2t(self.end)]
        
        return tuple(info)

class Actuation(object):
    def __init__(self, axis, id, stage):
        
        self.id = id
        self.start = rs.CurveStartPoint(axis)
        self.end = rs.CurveEndPoint(axis)
        self.stage = stage
        self.amp = 1
    
    def axis(self):
        
        return rs.AddLine(self.start, self.end)
    
    @staticmethod
    def model(typeName, axis, id, stage):
        
        converted = None
        
        for subType in Actuation.__subclasses__():
            if(typeName == subType.name):
               converted = subType(axis, id, stage)
               break
        
        return converted
    
    def outputInfo(self):
        
        info = []
        info += [self.name]
        info += [self.id]
        info += [p2t(self.start)]
        info += [p2t(self.end)]
        info += [self.amp]
        info += [self.stage.id]
        
        return tuple(info)

class ConstraintSpace(object):
    def __init__(self, info):
        self.axes = [rs.CreateVector(vec) for vec in info[1]]
        self.refPt = rs.CreateVector(info[2])
        self.spanVecs = [rs.CreateVector(vec) for vec in info[3]]
        self.spanDeg = len(self.spanVecs)
        self.axesDeg = len(self.axes)
        
    def sample(self, spanCount=5, stepSize=5, rotCount=8):
        
        # sample anchor point
        anchors = self.sampleAnchors(spanCount, stepSize)
        # smaple vector
        vectors = self.sampleVectors(rotCount)
        
        return vectors, anchors
        
    def sampleAnchors(self, spanCount, stepSize):
        
        stepCount = spanCount - 1
        
        if(self.spanDeg == 0):
            return self.refPt
        
        anchors = []
        for i in range(stepCount):
            xSteps = i - stepCount / 2
            xVec = rs.VectorScale(self.spanVecs[0], xSteps * stepSize)
            if(self.spanDeg == 1):
                anchor = rs.VectorAdd(self.refPt, xVec)
                anchors += [anchor]
            else:
                for j in range(stepCount): 
                    ySteps = j - stepCount / 2
                    yVec = rs.VectorScale(self.spanVecs[1], ySteps * stepSize)
                    if(self.spanDeg == 2):
                        vec = rs.VectorAdd(xVec, yVec)
                        anchor = rs.VectorAdd(self.refPt, vec)
                        anchors += [anchor]
                    else:
                        for k in range(stepCount): 
                            zSteps = k - stepCount / 2
                            zVec = rs.VectorScale(self.spanVecs[1], zSteps * stepSize)
                            vec = rs.VectorAdd(rs.VectorAdd(xVec, yVec), zVec)
                            anchor = rs.VectorAdd(self.refPt, vec)
                            anchors += [anchor]
        
        return anchors
    
    def sampleVectors(self, rotCount):
        
        unitStep = 360 / rotCount
        
        vectors = []
        if(self.axesDeg == 1):
            vec = self.axes[0]
            vectors += [vec]
        else:
            for i in range(rotCount):
                xFac = math.cos(math.radians(i * unitStep))
                vec = rs.VectorScale(self.axes[0], xFac)
                for j in range(rotCount):
                    yFac = math.cos(math.radians(j * unitStep))
                    yVec = rs.VectorScale(self.axes[1], yFac)
                    vec = rs.VectorAdd(vec, yVec)
                    if(self.axesDeg == 2):
                        vectors += [vec]
                    else:
                        for k in range(rotCount):
                            zFac = math.cos(math.radians(k * unitStep))
                            zVec = rs.VectorScale(self.axes[2], zFac)
                            vec = rs.VectorAdd(vec, zVec)
                            vectors += [vec]
        
        result = []
        for vec in vectors:
            if(rs.VectorLength(vec) < EPSILON):
                continue
            result += [rs.VectorUnitize(vec)]
        
        return result

################################################################################
# functional classes
################################################################################

class Graph(object):
    def __init__(self, nodeCount, edgeCount, adjacency):
        self.nodeCount = nodeCount
        self.edgeCount = edgeCount
        self.adj = self.adjacency = adjacency
        
    def __repr__(self):
        
        msg = "<Graph #N%d #E%d" % (self.nodeCount, self.edgeCount)
        if(DEBUG_MODE):
            msg += " adjs: %s" % str(self.adj)
        
        msg += '>'
        
        return msg
    
    def outputInfo(self):
        
        info = []
        info += [self.nodeCount]
        info += [self.edgeCount]
        info += [tuple(self.adj)]
        
        return tuple(info)

class Modeler(object):
    def __init__(self):
        self.cons = []
        self.stages = []
        self.goals = []
        self.acts = []
    
    def reset(self):
        self.cons = []
        self.stages = []
        self.goals = []
        self.acts = []
    
    def addWires(self):
        # adds a wire element
        
        # select wire
        prompt = "Select lines representing wires"
        filter = Rhino.DocObjects.ObjectType.Curve
        wires = Modeler.getObjects(prompt, filter)
        lines = []
        isValid = False
        for result in wires:
            if(result != None):
                line = result.ObjectId
                if(line not in self.cons):
                    lines += [line]
                    isValid = True
        
        # select first stage
        prompt = "select first connected stage"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage0 = Modeler.getObjects(prompt, filter, True)
        s0 = stage0.ObjectId
        
        # select second stage
        prompt = "select second connected stage"
        stage1 = Modeler.getObjects(prompt, filter, True)
        s1 = stage1.ObjectId
        
        if(len(lines) > 0 and isValid): # valid input exists
            # add stages to contrainer (no duplicates)
            if(s0 not in self.stages):
                self.stages += [InputObject(Stage.name, s0)]
            if(s1 not in self.stages):
                self.stages += [InputObject(Stage.name, s1)]
            
            # create constraint elements (no duplicates)
            s0Index = self.stages.index(s0)
            s1Index = self.stages.index(s1)
            adjs = sorted([s0Index, s1Index])
            adjs = tuple([self.stages[adjs[0]], self.stages[adjs[1]]]) 
            
            self.cons += [InputObject(Wire.name, line, adjs) for line in lines]
    
    def addPlaceHolder(self):
        
        # select first stage
        prompt = "select first connected stage"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage0 = Modeler.getObjects(prompt, filter, True)
        s0 = stage0.ObjectId
        
        # select second stage
        prompt = "select second connected stage"
        stage1 = Modeler.getObjects(prompt, filter, True)
        s1 = stage1.ObjectId
        
        if(s0 and s1):
            # add stages to contrainer (no duplicates)
            if(s0 not in self.stages):
                self.stages += [InputObject(Stage.name, s0)]
            if(s1 not in self.stages):
                self.stages += [InputObject(Stage.name, s1)]
                
            # create constraint elements (no duplicates)
            s0Index = self.stages.index(s0)
            s1Index = self.stages.index(s1)
            adjs = sorted([s0Index, s1Index])
            adjs = tuple([self.stages[adjs[0]], self.stages[adjs[1]]]) 
            
            if(adjs not in self.cons):
                self.cons += [InputObject(PlaceHolder.name, adj=adjs)]
    
    def addGoalTrans(self):
        
        # select axes
        prompt = "Select lines representing translations"
        filter = Rhino.DocObjects.ObjectType.Curve
        axes = Modeler.getObjects(prompt, filter)
        lines = []
        isValid = False
        for result in axes:
            if(result != None): # valid input check
                axis = result.ObjectId
                if(InputObject(GoalTrans.name, axis) not in self.goals): # duplicate check
                    lines += [axis]
                    isValid = True
        
        # select stage
        prompt = "select the stage these DOFs are attached to"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage = Modeler.getObjects(prompt, filter, True)
        s = stage.ObjectId
        
        if(len(lines) > 0 and isValid): # valid input exists
            # add stages to contrainer (no duplicates)
            if(s not in self.stages):
                self.stages += [InputObject(Stage.name, s)]
            
            # create goal DOFs (no duplicates)
            sIndex = self.stages.index(s)
            adjs = self.stages[sIndex]
            self.goals += [InputObject(GoalTrans.name, line, adjs) for line in lines]
        
    def addGoalRot(self):
        
        # select axes
        prompt = "Select lines representing rotation axes"
        filter = Rhino.DocObjects.ObjectType.Curve
        axes = Modeler.getObjects(prompt, filter)
        lines = []
        isValid = False
        for result in axes:
            if(result != None): # valid input check
                axis = result.ObjectId
                if(InputObject(GoalRots.name, axis) not in self.goals): # duplicate check
                    lines += [axis]
                    isValid = True
        
        # select stage
        prompt = "select the stage these DOFs are attached to"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage = Modeler.getObjects(prompt, filter, True)
        s = stage.ObjectId
        
        if(len(lines) > 0 and isValid): # valid input exists
            # add stages to contrainer (no duplicates)
            if(s not in self.stages):
                self.stages += [InputObject(Stage.name, s)]
            
            # create goal DOFs (no duplicates)
            sIndex = self.stages.index(s)
            adjs = self.stages[sIndex]
            self.goals += [InputObject(GoalRots.name, line, adjs) for line in lines]
    
    def addActTrans(self):
        
        # select axis
        prompt = "Select lines representing translations"
        filter = Rhino.DocObjects.ObjectType.Curve
        axis = Modeler.getObjects(prompt, filter, True)
        line = axis.ObjectId
        if(InputObject(ActTrans.name, line) not in self.acts):
            isValid = True
        
        # select stage
        prompt = "select the stage this actuation is attached to"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage = Modeler.getObjects(prompt, filter, True)
        s = stage.ObjectId
        
        if(isValid): # valid input exists
            # add stages to contrainer (no duplicates)
            if(s not in self.stages):
                self.stages += [InputObject(Stage.name, s)]
            
            # create goal DOFs (no duplicates)
            sIndex = self.stages.index(s)
            adjs = self.stages[sIndex]
            self.acts += [InputObject(ActTrans.name, line, adjs)]
    
    def addActRot(self):
        
        # select axis
        prompt = "Select line representing rotation axis"
        filter = Rhino.DocObjects.ObjectType.Curve
        axis = Modeler.getObjects(prompt, filter, True)
        line = axis.ObjectId
        if(InputObject(ActRot.name, line) not in self.acts):
            isValid = True
        
        # select stage
        prompt = "select the stage this actuation is attached to"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        stage = Modeler.getObjects(prompt, filter, True)
        s = stage.ObjectId
        
        if(isValid): # valid input exists
            # add stages to contrainer (no duplicates)
            if(s not in self.stages):
                self.stages += [InputObject(Stage.name, s)]
            
            # create goal DOFs (no duplicates)
            sIndex = self.stages.index(s)
            adjs = self.stages[sIndex]
            self.acts += [InputObject(ActRot.name, line, adjs)]
    
    def deleteObjects(self):
        
        # select stage
        prompt = "select objects to remove from the model"
        filter = Rhino.DocObjects.ObjectType.AnyObject
        objects = Modeler.getObjects(prompt, filter)
        objects = [obj.ObjectId for obj in objects]
        
        # flag objects to be deleted
        for obj in objects:
            if(obj in self.cons):
                self.cons[self.cons.index(obj)].deleteFlag = True
            if(obj in self.stages):
                self.stages[self.stages.index(obj)].deleteFlag = True
            if(obj in self.goals):
                self.goals[self.goals.index(obj)].deleteFlag = True
            if(obj in self.acts):
                self.acts[self.acts.index(obj)].deleteFlag = True
    
    def checkObjects(self):
        
        delTypes = [0, 0, 0]
        
        validCons = []
        for con in self.cons:
            isValid, delType = con.checkValidity()
            if(not isValid):
                delTypes[delType] += 1
            else:
                validCons += [con]
        
        validStages = []
        for stage in self.stages:
            isValid, delType = stage.checkValidity()
            if(not isValid):
                delTypes[delType] += 1
            else:
                validStages += [stage]
            
        validGoals = []
        for goal in self.goals:
            isValid, delType = goal.checkValidity()
            if(not isValid):
                delTypes[delType] += 1
            else:
                validGoals += [goal]
        
        self.cons = validCons
        self.stages = validStages
        self.goals = validGoals
        
        msg = ''
        if(delTypes[0] != 0):
            msg += "%d object(s) were deleted by user commands. " % delTypes[0]
        if(delTypes[1] != 0):
            msg += "%d object(s) were deleted because their reference objects were removed. " % delTypes[1]
        if(delTypes[2] != 0):
            msg += "%d object(s) were deleted becasue their associated objects were removed. " % delTypes[2]
        if(len(msg) == 0):
            msg = "Model was successfully processed."
        
        return msg

    @staticmethod
    def getObjects(prompt, filter, single=False):
    
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        if(single):
            result = Rhino.Input.RhinoGet.GetOneObject(prompt, False, filter)
        else:
            result = Rhino.Input.RhinoGet.GetMultipleObjects(prompt, False, filter)
        sc.doc = ghdoc
        
        return result[1]

class ModelerSimp(object):
    def __init__(self):
        self.cons = []
        self.baseStage = None
        self.freeStage = None
        self.goals = []
    
    def resetStages(self):
        
        self.baseStage = None
        self.freeStage = None
        
    def addStage(self, targ):

        # get object
        prompt = "Select a solid representing the %s stage." % targ
        stage, success = self._addObject(prompt)
        
        # generate message based on input and condition
        if(not success):
            msg = "Something went wrong. Please try again."
        else:
            stage = InputObject(Stage.name, stage)
            msg = "Successful!"
            if(targ == "base"):
                self.baseStage = stage
            elif(targ == "free"):
                self.freeStage = stage

        # stages status check
        if(self.baseStage != None and self.freeStage != None):
            msg += " All stages are set."
        elif(self.baseStage != None and self.freeStage == None):
            msg += " Still need to set the free stage."
        elif(self.baseStage == None and self.freeStage != None):
            msg += " Still need to set the base stage."
        else:
            msg += " Both stages are not specified."
        
        return msg

    def addGoal(self, mode):

        if(mode == GoalTrans.nameShort):
            promptName = GoalTrans.promptName
            tag = GoalTrans.name
        elif(mode == GoalRots.nameShort):
            promptName = GoalRots.promptName
            tag = GoalRots.name

        lines, success = self._addLines(promptName)

        if(success):
            goals = [InputObject(tag, line) for line in lines]
            self.goals += goals

    def resetCons(self):

        self.cons = []

    def addDelRod(self, delete=False):
        
        self._addDelCons(Wire.promptName, Wire.name, delete)

    def addDelCable(self, delete=False):

        self._addDelCons(Cable.promptName, Cable.name, delete)

    def addDelSensor(self, delete=False):

        self._addDelCons(Sensor.promptName, Sensor.name, delete)

    def printInfo(self, mode=None):

        if(mode == None):
            print(self.cons, self.baseStage, self.freeStage, self.goals)
        elif(mode == "shared"):
            print(self.cons, self.baseStage, self.freeStage)
        elif(mode == "goals"):
            print(self.goals)
        elif(mode == "cables"):
            print(self.cons)
        elif(mode == "sensors"):
            print(self.cons)

    def _addDelCons(self, promptName, tag, delete=False):

        if(delete): promptName += " to delete"
        
        lines, success = self._addLines(promptName)
        
        if(success):
            cons = [InputObject(tag, line) for line in lines]
            # avoiding duplicates
            for con in cons:
                if(not delete and con not in self.cons):
                    self.cons += [con]
                elif(delete and con in self.cons):
                    self.cons.pop(self.cons.index(con))

    def _addLines(self, lineName="wires"):
        # adds a wire element
        
        # prompt
        prompt = "Select lines representing %s" % lineName
        filter = Rhino.DocObjects.ObjectType.Curve
        wires = Modeler.getObjects(prompt, filter)
        
        # clean
        lines = []
        isValid = False
        for result in wires:
            if(result != None):
                lines += [result.ObjectId]
                isValid = True
        
        return lines, isValid
        
    def _addObjects(self, prompt="Select objects to add"):
        # adds a stage object
        filter = Rhino.DocObjects.ObjectType.AnyObject
        selObjects = Modeler.getObjects(prompt, filter)
        
        # clean
        objects = []
        isValid = False
        for obj in selObjects:
            if(obj != None):
                objects += obj.ObjectId
                isValid = True

        return objects, isValid
    
    def _addObject(self, prompt="Select object to add"):
        # adds a stage object
        filter = Rhino.DocObjects.ObjectType.AnyObject
        selObject = Modeler.getObjects(prompt, filter, True)
        
        if(selObject.ObjectId != None):
            isValid = True
            selObject = selObject.ObjectId
        else:
            isValid = False
            selObject = None
            
        return selObject, isValid

    @staticmethod
    def getObjects(prompt, filter, single=False):
    
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        if(single):
            result = Rhino.Input.RhinoGet.GetOneObject(prompt, False, filter)
        else:
            result = Rhino.Input.RhinoGet.GetMultipleObjects(prompt, False, filter)
        sc.doc = ghdoc
        
        return result[1]

class Model(object):
    def __init__(self, modeler, mats):
        self.c = self.constraints = []
        self.s = self.stages = []
        self.goals = []
        self.acts = []
        self.g = self.graph = None
        self.m = self.materials = mats
        self.modeler = modeler
        
        self._model()
    
    def _model(self):
        # init function
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        self._modelStages()
        self._modelConstraints()
        self._modelJoints()
        self._modelGraph()
        self._modelGoals()
        self._modelActs()
        sc.doc = ghdoc
        if(DEBUG_MODE):
            self.printInfo()
        
    def _modelStages(self):
        
        stages = []
        for i in range(len(self.modeler.stages)):
            s = self.modeler.stages[i]
            stage = Stage.model(s.type, s.geo, i)
            assert stage != None, MSG_UNKNOWN_STAGE_TYPE
            stages += [stage]
        
        self.s = self.stages = stages
    
    def _modelConstraints(self):
        
        # real constraint elements
        constraints = []
        for i in range(len(self.modeler.cons)):
            c = self.modeler.cons[i]
            adjSIds = [self.modeler.stages.index(s) for s in c.adj]
            adjStages = [self.s[id] for id in adjSIds]
            elem = Constraint.model(c.type, c.geo, i, adjStages)
            assert elem != None, MSG_UNKNOWN_CON_TYPE
            constraints += [elem]
        
        self.c = self.constraints = constraints
        
    def _modelJoints(self):
        
        # joints from constraints
        jMap = {}
        for con in self.c:
            adj = tuple([s.id for s in con.adjs])
            jMap[adj] = jMap.get(adj, []) + [con]
        
        joints = []
        adjs = jMap.keys()
        for i in range(len(adjs)):
            adj = adjs[i]
            cons = jMap[adj]
            stages = [self.s[sId] for sId in adj]
            joint = Joint(cons, i, stages)
            joints += [joint]
        
        self.j = self.joints = joints
       
    def _modelGraph(self):
        
        conMap = {}
        for c in self.c:
            conMap[c.id] = [s.id for s in c.adjs]
        
        adjs = []
        for id in sorted(conMap.keys()):
            adjs += [conMap[id]]
        
        graph = Graph(len(self.s), len(self.c), adjs)
        self.g = self.graph = graph
    
    def _modelGoals(self):
        
        goals = []
        for i in range(len(self.modeler.goals)):
            g = self.modeler.goals[i]
            stage = self.s[self.modeler.stages.index(g.adj)]
            axis = g.geo
            goal = Goal.model(g.type, axis, i, stage)
            goals += [goal]
        
        self.goals = goals
    
    def _modelActs(self):
        
        acts = []
        for i in range(len(self.modeler.acts)):
            act = self.modeler.acts[i]
            stage = self.s[self.modeler.stages.index(act.adj)]
            axis = act.geo
            actuation = Actuation.model(act.type, axis, i, stage)
            acts += [actuation]
        
        self.acts = acts
    
    def setFixedStage(self, id, isFixed):
        
        if(id >= len(self.s)):
            return
        
        for s in self.s:
            s.isFixed = False if s.id != id else True
        
        self.s[id].isFixed=isFixed
    
    def printInfo(self):
        
        print("Printing model info")
        print("Constaints:") 
        for c in self.c:
            print(c)
        print("Joints:") 
        for j in self.j:
            print(j)
        print("Stages:")
        for s in self.s:
            print(s)
        print("Graph:")
        print(self.g)
    
    def placeHolderPt(self):
        # shortcut to the first one found
        
        for c in self.c:
            if(isinstance(c, PlaceHolder)):
                return c.center
    
    def outputInfo(self):
        
        output = []
        output += [self._outputGraphInfo()]
        output += [self._outputMatsInfo()]
        output += [self._outputConsInfo()]
        output += [self._outputJointsInfo()]
        output += [self._outputStagesInfo()]
        output += [self._outputGoalsInfo()]
        output += [self._outputActsInfo()]
        
        return tuple(output)
    
    def _outputGraphInfo(self):
        
        return self.g.outputInfo()
    
    def _outputMatsInfo(self):
        
        output = []
        for m in self.m:
            output += [m.outputInfo()]
        
        return tuple(output)
    
    def _outputConsInfo(self):
        
        output = []
        for c in self.c:
            output += [c.outputInfo()]
        
        return tuple(output)
    
    def _outputJointsInfo(self):
        
        output = []
        for j in self.j:
            output += [j.outputInfo()]
        
        return tuple(output)
    
    def _outputStagesInfo(self):
        
        output = []
        for s in self.s:
            output += [s.outputInfo()]
        
        return tuple(output)
    
    def _outputGoalsInfo(self):
        
        output = []
        for g in self.goals:
            output += [g.outputInfo()]
        
        return tuple(output)
    
    def _outputActsInfo(self):
        
        output = []
        for act in self.acts:
            output += [act.outputInfo()]
        
        return tuple(output)
        
    def outputTags(self, palette, viewCons, viewJoint, viewStage):
        
        pts, ids, clrs = [], [], []
        
        if(viewCons):
            for c in self.c:
                pt, id = c.outputTags()
                clr = palette['c']
                pts += [pt]
                ids += [id]
                clrs += [clr]
        
        if(viewJoint):
            for j in self.j:
                pt, id = j.outputTags()
                clr = palette['j']
                pts += [pt]
                ids += [id]
                clrs += [clr]
        
        if(viewStage):
            for s in self.s:
                pt, id = s.outputTags()
                clr = palette['s']
                pts += [pt]
                ids += [id]
                clrs += [clr]
        
        return pts, ids, clrs

class ModelSimp(object):
    def __init__(self, CMCore):

        self.core = CMCore

        self.modes = []
        self.stages = []
        self.cables = []
        self.sensors = []
        self.rods = []
        
        self._model()

    def _model(self):
        
        # init function
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        self._modelStages()
        self._modelModes()
        self._modelRods()
        self._modelCables()
        self._modelSensors()
        sc.doc = ghdoc

    def _modelStages(self):

        modeler = self.core.modelerShared
        source = []
        if(modeler.baseStage != None):
            source += [modeler.baseStage]
        if(modeler.freeStage != None):
            source += [modeler.freeStage]

        stages = []
        for i in range(len(source)):
            s = source[i]
            stage = Stage.model(s.type, s.geo, i)
            assert stage != None, MSG_UNKNOWN_STAGE_TYPE
            stages += [stage]
        
        for i in range(2 - len(source)):
            stage = Stage(None, i + len(source), len(stages) == 0)
            stages += [stage]

        self.stages = stages

    def _modelModes(self):
        
        modelers = self.core.modelerGoals
        targ = 1

        goals = []
        for i in range(len(modelers)):
            temp = []
            for j in range(len(modelers[i].goals)):
                g = modelers[i].goals[j]
                adjStage = self.stages[targ]
                goal = Goal.model(g.type, g.geo, j, adjStage)
                assert goal != None, MSG_UNKNOWN_CON_TYPE
                temp += [goal]
            goals += [temp]
        
        self.modes = goals

    def _modelRods(self):
        
        modeler = self.core.modelerShared
        adjs = (0, 1)

        # real constraint elements
        constraints = []
        for i in range(len(modeler.cons)):
            c = modeler.cons[i]
            adjStages = [self.stages[id] for id in adjs]
            elem = Constraint.model(c.type, c.geo, i, adjStages)
            assert elem != None, MSG_UNKNOWN_CON_TYPE
            constraints += [elem]
        
        self.rods = constraints

    def _modelCables(self):
        
        modelers = self.core.modelerCables
        adjs = (0, 1)

        # real constraint elements
        constraints = []
        for i in range(len(modelers)):
            temp = []
            for j in range(len(modelers[i].cons)):
                c = modelers[i].cons[j]
                adjStages = [self.stages[id] for id in adjs]
                elem = Constraint.model(c.type, c.geo, j, adjStages)
                assert elem != None, MSG_UNKNOWN_CON_TYPE
                temp += [elem]
            constraints += [temp]
        
        self.cables = constraints

    def _modelSensors(self):
        
        modelers = self.core.modelerSensors
        adjs = (0, 1)

        # real constraint elements
        constraints = []
        for i in range(len(modelers)):
            temp = []
            for j in range(len(modelers[i].cons)):
                c = modelers[i].cons[j]
                adjStages = [self.stages[id] for id in adjs]
                elem = Constraint.model(c.type, c.geo, j, adjStages)
                assert elem != None, MSG_UNKNOWN_CON_TYPE
                temp += [elem]
            constraints += [temp]
        
        self.sensors = constraints

    def outputInfo(self):
        
        output = []
        output += [self._outputStagesInfo()]
        output += [self._outputModesInfo()]
        output += [self._outputRodsInfo()]
        output += [self._outputCablesInfo()]
        output += [self._outputSensorsInfo()]
        
        return tuple(output)

    def _outputStagesInfo(self):
        
        output = []
        for s in self.stages:
            output += [s.outputInfo()]
        
        return tuple(output)

    def _outputModesInfo(self):
        
        output = []
        for mode in self.modes:
            temp = []
            for goal in mode:
                temp += [goal.outputInfo()]
            output += [tuple(temp)]
        
        return tuple(output)

    def _outputRodsInfo(self):
        
        output = []
        for rod in self.rods:
            output += [rod.outputInfo()]
        
        return tuple(output)
    
    def _outputCablesInfo(self):

        output = []
        for group in self.cables:
            temp = []
            for cable in group:
                temp += [cable.outputInfo()]
            output += [tuple(temp)]
        
        return tuple(output)
    
    def _outputSensorsInfo(self):
        
        output = []
        for group in self.sensors:
            temp = []
            for sensor in group:
                temp += [sensor.outputInfo()]
            output += [tuple(temp)]
        
        return tuple(output)
    
################################################################################
# specific constraint types
################################################################################

class Wire(Constraint):
    name = "WireCircular"
    promptName = "flexural rods"
    degree = 5
    
    def __init__(self, line, id, adjs):
        
        super(Wire, self).__init__(line, id, adjs)
        
        self.start = rs.CurveStartPoint(line)
        self.end = rs.CurveEndPoint(line)
        self.len = rs.Distance(self.start, self.end)
        self.vec = rs.VectorUnitize(rs.VectorSubtract(self.end, self.start))
        self.rad = self.len / 40 # hard-coded rod radius
        
    def __repr__(self):
        
        msg = "<Wire #%d" % self.id
        if(DEBUG_MODE):
            msg += " adj:%s s:(%s) e:(%s)" % (str(self.adjs), str(self.start), str(self.end))
            
        msg += '>'
        
        return msg
    
    def _tagPt(self):
        
        pt = rs.PointAdd(self.start, self.end)
        pt = rs.PointDivide(pt, 2)
        
        return pt
    
    def outputInfo(self):
        
        info = super(Wire, self).outputBasicInfo()
        info += [p2t(self.start)]
        info += [p2t(self.end)]
        info += [str(dFix(self.rad))]
        
        return tuple(info)
    
    def outputVis(self):
        
        return rs.AddPipe(rs.AddLine(self.start, self.end), 0, self.rad, cap=2)
    
    def frames(self):
        
        p0 = rs.PlaneFromNormal(self.start, self.vec)
        p1 = rs.PlaneFromNormal(self.end, rs.VectorReverse(self.vec))
        
        return p0, p1
    
    def findNeighborStage(self):
        
        s0, s1 = self.adjs
        startDistS0 = s0.minDistTo(self.start)
        startDistS1 = s1.minDistTo(self.start)
        
        if(startDistS0 > startDistS1):
            return s1.id, s0.id
        else:
            return s0.id, s1.id
    
    @staticmethod
    def _isType(geometry):
        
        isLine = rs.IsLine(geometry)
        
        return isLine

class Sensor(Wire):
    name = "Sensor"
    promptName = "stretchable sensors"

    def __init__(self, line, id, adjs):

        super(Sensor, self).__init__(line, id, adjs)
        self.rad = 1.5 # hard-coded sensor radius

class Cable(Wire):
    name = "Cable"
    promptName = "tensioning cables"

    def __init__(self, line, id, adjs):

        super(Cable, self).__init__(line, id, adjs)

        self.rad = .25 # hard-coded cable radius
    
class PlaceHolder(Constraint):
    name = "PlaceHolder"
    degree = 6
    
    def __init__(self, geo, id, adjs):
        
        super(PlaceHolder, self).__init__(geo, id, adjs)
        self.center = None
        
        p0 = self.adjs[0].center
        p1 = self.adjs[1].center
        self.center = rs.PointDivide(rs.PointAdd(p0, p1), 2)
        self.conSpaces = []
        self.consCandidate = []
    
    def __repr__(self):
        
        msg = "<PlaceHolder #%d" % self.id
        if(DEBUG_MODE):
            msg += " adj:%s s:(%s) e:(%s)" % (str(self.adjs), str(self.start), str(self.end))
            
        msg += '>'
        
        return msg
    
    def outputInfo(self):
        
        info = super(PlaceHolder, self).outputBasicInfo()
        info += [p2t(self.center)]
        
        return tuple(info)
    
    def outputVis(self):
        return [rs.AddSphere(self.center, DEFAULT_WIRE_RAD)]
    
    def _tagPt(self):
        
        return self.center
    
    def addConSpaces(self, conSpaces):
        
        for conSpace in conSpaces:
            self.conSpaces += [ConstraintSpace(conSpace)]
        
    def sampleConstraint(self, id=0):
        
        self.consCandidate = []
        
        vecs, pts = self.conSpaces[id].sample()
        valids, invalids = [], []
        
        
        for vec in vecs:
            for pt in pts:
                isValid = self.checkConstraint(pt, vec)
                if(not isValid):
                    invalids += [PlaceHolder._addLine(pt, vec, 3)]
        
        for pts in self.consCandidate:
            valids += [rs.AddLine(pts[0], pts[1])]
        
        return valids, invalids
    
    def checkConstraint(self, vec, pt, length=1e+6):
        
        sc.doc = Rhino.RhinoDoc.ActiveDoc
        
        line = PlaceHolder._addLine(vec, pt)
        int0 = self.adjs[0].intersect(line)
        int1 = self.adjs[1].intersect(line)
        if(int0 != None and int1 != None):
            isValid = True
            self.addConsCandidate([int0, int1])
        else:
            isValid = False
        
        rs.DeleteObject(line)
        
        sc.doc = ghdoc
        
        return isValid
    
    def addConsCandidate(self, newPts):
        
        sc.doc = ghdoc
        
        isNovel = True
        for pts in self.consCandidate:
            dist0 = rs.Distance(pts[0], newPts[0])
            dist1 = rs.Distance(pts[1], newPts[1])
            if(dist0 < EPSILON and dist1 < EPSILON):
                isNovel = False
                break
        
        if(isNovel):
            self.consCandidate += [newPts]
                
        sc.doc = Rhino.RhinoDoc.ActiveDoc
    
    @staticmethod
    def _addLine(pt, vec, length=1e+6):
        
        start = rs.VectorAdd(pt, rs.VectorScale(vec, -.5 * length))
        end = rs.VectorAdd(pt, rs.VectorScale(vec, .5 * length))
        line = rs.AddLine(start, end)
        
        return line
        
    @staticmethod
    def _isType(geometry):
        
        return geometry == None

class GoalTrans(Goal):
    name = "GoalTrans"
    nameShort = "trans"
    promptName = "translation axis"
    
    def __init__(self, *args):
        
        if(len(args) == 3):
            axis, id, stageId = args
        else:
            axis, id = args
            stageId = 1

        super(GoalTrans, self).__init__(axis, id, stageId)
        self.type = 't'
    
    def outputVis(self):
        
        return drawTransAxis(self.axis())

class GoalRots(Goal):
    name = "GoalRots"
    nameShort = "rot"
    promptName = "rotation axis"
    
    def __init__(self, *args):
        print(args)
        if(len(args) == 3):
            axis, id, stageId = args
        else:
            axis, id = args
            stageId = 1

        super(GoalRots, self).__init__(axis, id, stageId)
        self.type = 'r'
    
    def outputVis(self):
        
        return drawRotAxis(self.axis())

class ActTrans(Actuation):
    name = "ActTrans"
    
    def __init__(self, axis, id, stage):
        super(ActTrans, self).__init__(axis, id, stage)
        self.type = 't'

class ActRot(Actuation):
    name = "ActRot"
    
    def __init__(self, axis, id, stage):
        super(ActRot, self).__init__(axis, id, stage)
        self.type = 'r'

################################################################################
# utilities
################################################################################

def p2t(point):
    # convert a point coordinate into a tuple
        
    return tuple([str(dFix(num)) for num in point])

def dFix(num):
    # decimal rounding
    
    return Decimal(num) + Decimal(0)

def drawArrow(anchor, vec):
    
    scale = rs.VectorLength(vec) / SCALE
    
    height = 2 * scale
    dia = 1 * scale
    
    end = rs.VectorAdd(anchor, vec)
    pipe = rs.AddPipe(rs.AddLine(anchor, end), 0, dia / 4, cap=1)
    apex = rs.VectorAdd(end, rs.VectorScale(rs.VectorUnitize(vec), height))
    plane = rs.PlaneFromNormal(apex, -vec)
    
    cone = rs.AddCone(plane, height, dia)
    arrow = rs.BooleanUnion(pipe + [cone])
    
    return arrow[0]

def drawTransAxis(line):
    
    start = rs.CurveStartPoint(line) 
    mid = rs.CurveMidPoint(line)
    end = rs.CurveEndPoint(line)
    vecMS = rs.PointSubtract(start, mid)
    vecME = rs.PointSubtract(end, mid)
    
    arrows = [drawArrow(mid, vecMS), drawArrow(mid, vecME)]
    
    return arrows

def drawRotAxis(line, drawArc=True):
    
    scale = rs.CurveLength(line) / SCALE
    
    height = 2 * scale
    dia = 1 * scale
    arcDeg = 240
    rotation = (360 - arcDeg) / 2
    arcRad = 4 * scale
    
    start = rs.CurveStartPoint(line) 
    mid = rs.CurveMidPoint(line)
    end = rs.CurveEndPoint(line)
    vec = rs.VectorUnitize(rs.VectorSubtract(end, start))
    
    axisPipe = rs.AddPipe(line, 0, dia / 4, cap=2)
    arcPlane = rs.RotatePlane(rs.PlaneFromNormal(mid, vec), rotation, vec)
    arc = rs.AddArc(arcPlane, arcRad, arcDeg)
    arcPipe = rs.AddPipe(arc, 0, dia / 4, cap=1)
    
    arcDom = rs.CurveDomain(arc)
    tangStart = rs.CurveTangent(arc, arcDom[0])
    tangEnd   = rs.CurveTangent(arc, arcDom[1]) * -1
    arcStart = rs.CurveStartPoint(arc)
    arcEnd = rs.CurveEndPoint(arc)
    
    planeStart = rs.PlaneFromNormal(arcStart - tangStart, tangStart)
    planeEnd = rs.PlaneFromNormal(arcEnd - tangEnd, tangEnd)
    coneStart = rs.AddCone(planeStart, height, dia)
    coneEnd = rs.AddCone(planeEnd, height, dia)
    arcArrow = rs.BooleanUnion(arcPipe + [coneStart, coneEnd])
    
    if(drawArc): shape = arcArrow + axisPipe
    else: shape = axisPipe
    
    return shape

def exhaustiveSubclasses(targClass):
    result = [targClass]

    subclasses = targClass.__subclasses__()
    if(len(subclasses) != 0):
        for subclass in subclasses:
            moreSubclassses = exhaustiveSubclasses(subclass)
            result += moreSubclassses

    return result

################################################################################
# main
################################################################################

class CMCore(object):
    def __init__(self, simpMode=False):
        self.modeler = Modeler()
        self.model = None
        self.materials = []
        self.tarId = None
        self.groundId = None
        
        if(simpMode):
            self.modelerShared = ModelerSimp()
            self.modelerGoals = [ModelerSimp()]
            self.modelerCables = []
            self.modelerSensors = [ModelerSimp()]
    
        self.modeCount = len(self.modelerGoals)
        self.sensorCount = len(self.modelerSensors)
        self.cableGroupCount = 0
        self.conSubspaceCount = 0
        self.modeFreeSpaceCount = '0'
        self.cableConSpaceCount = '0'
        self.sensorConSpaceCount = '0'
        
        self.snapShot = None

        self.msgStages = "Both stages are not specified."
        
    def createFromModeler(self):
        
        self.model = Model(self.modeler, self.materials)
    
    def createMaterial(self, name, E, G, D, MS):
        
        id = len(self.materials)
        material = IsoMat(name, id, E, G, D, MS)
        if(material not in self.materials):
            self.materials += [material]
    
    def outputInfo(self):
        
        if(isinstance(self.model, Model) or isinstance(self.model, ModelSimp)):
            return str(self.model.outputInfo())
    
    def outputTags(self, palette, viewConsId, viewJointId, viewStageId):
        
        if(isinstance(self.model, Model)):
            return self.model.outputTags(palette, viewConsId, viewJointId, viewStageId)
    
    def outputGeos(self):
        cons = [c.geo for c in self.modeler.cons]
        stages = [s.geo for s in self.modeler.stages]
        
        return cons, stages
    
    def outputViewCenter(self, mode, stageId=None):
        
        if(not isinstance(self.model, Model)): return # uninitialized shortcut
        
        if(mode == 0): # target stage
            if(self.tarId < len(self.model.s)):
                pt, _ = self.model.s[self.tarId].outputTags()
                return pt
        elif(mode == 1): # placeholder
            return self.model.placeHolderPt()
             
        else:
            assert False, MSG_UNIMPLEMENTED
    
    def outputVis(self, args):
        
        if(not isinstance(self.model, Model)): # uninitialized shortcut
            return [None] * len(args)
        
        result = []
        for arg in args:
            shapes = []
            
            if(arg == 't'): # target stage
                if(self.tarId < len(self.model.s)):
                    shapes += self.model.s[self.tarId].outputVis()
                
            elif(arg == 'g'): # ground stage
                if(self.groundId < len(self.model.s)):
                    shapes += self.model.s[self.groundId].outputVis()
                
            elif(arg == 'o'): # other stage
                for i in range(len(self.model.s)):
                    if(i == self.groundId or i == self.tarId): continue
                    shapes += self.model.s[i].outputVis()
                
            elif(arg == 'c'): # constraints
                for c in self.model.c:
                    shapes += c.outputVis()
                
            else:
                assert False, MSG_UNIMPLEMENTED
            
            result += [shapes]
        
        return result
    
    def outputGoals(self):
        
        trans = []
        rots = []
        for goal in self.model.goals:
            
            if(isinstance(goal, GoalTrans)):
                trans += goal.outputVis()
            else:
                rots += goal.outputVis()
        
        return rots, trans
        
    def setFixedStage(self, id, isFixed=True):
        
        self.model.setFixedStage(id, isFixed)

    def test(self, motions, id):
        self.model.c[0].addConSpaces(motions)
        valids, invalids = self.model.c[0].sampleConstraint(id)
        
        return valids, invalids

    ################################################################################
    # simple mode functions

    def addBaseStage(self):
        self.msgStages = self.modelerShared.addStage("base")
    
    def addFreeStage(self):
        self.msgStages = self.modelerShared.addStage("free")

    def resetModes(self):
        
        self.modelerShared.resetStages()
        self.modelerGoals = [ModelerSimp()]
        self.modeCount = len(self.modelerGoals)

    def addMode(self):
        
        self.modelerGoals += [ModelerSimp()]
        self.modeCount = len(self.modelerGoals)
        
        self.modelerSensors += [ModelerSimp()]
        self.sensorCount = len(self.modelerSensors)
    
    def delMode(self, id):
        
        self.modelerGoals.pop(id)
        self.modeCount = len(self.modelerGoals)
        
        self.modelerSensors.pop(id)
        self.sensorCount = len(self.modelerSensors)
    
    def addGoalRot(self, id):

        self.modelerGoals[id].addGoal(GoalRots.nameShort)

    def addGoalTrans(self, id):
        
        self.modelerGoals[id].addGoal(GoalTrans.nameShort)
    
    def resetRod(self):

        self.modelerShared.resetCons()

    def addRod(self):

        self.modelerShared.addDelRod()

    def delRod(self):
        
        self.modelerShared.addDelRod(True)

    def addCab(self, id):

        self.modelerCables[id].addDelCable()

    def delCab(self, id):
        
        self.modelerCables[id].addDelCable(True)

    def addSen(self, id):
        
        self.modelerSensors[id].addDelSensor()

    def delSen(self, id):
        
        self.modelerSensors[id].addDelSensor(True)

    def printInfo(self):

        print("Shared")
        self.modelerShared.printInfo("shared")
        print("Goal DOF Modes")
        for m in self.modelerGoals: m.printInfo("goals")
        print("Cable Groups")
        for m in self.modelerCables: m.printInfo("cables")
        print("Sensors groups")
        for m in self.modelerSensors: m.printInfo("sensors")
    
    def createModel(self):

        self.model = ModelSimp(self)

    def setCableGroupCount(self, count):
        
        curCount = len(self.modelerCables)
        if(curCount < count): # need more
            for i in range(count - curCount):
                self.modelerCables += [ModelerSimp()]
        elif(curCount > count): # toom many
            self.modelerCables = self.modelerCables[:count]
        
        self.cableGroupCount = len(self.modelerCables)

CM = CMCore(SimpMode)
if(reset):
    CM = CMCore(SimpMode)
else:
    CM = CM