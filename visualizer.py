"""Provides a scripting component.
    Inputs:
        x: The x script variable
        y: The y script variable
    Output:
        a: The a output variable"""

__author__ = "MML"
__version__ = "2022.01.03"

import rhinoscriptsyntax as rs
from Grasshopper import DataTree
from Grasshopper.Kernel.Data import GH_Path as treePath
from Rhino.Geometry import GeometryBase

SCALE = 30
EPSILON = 1e-6

def unpackBases(bases):
    if(len(bases) == 0):
        return []
    if(not isinstance(bases[0], tuple)):
        return rs.AddPoint(bases)
    else:
        result = []
        for item in bases:
            result += [unpackBases(item)]
        return result

def unpackSpans(vecs):
    
    spanVecs = []
    
    if(len(vecs) == 0): return spanVecs
    for vec in vecs:
        if(len(vec) == 0): continue
        newVec = rs.AddPoint(vec)
        spanVecs += [newVec]
    
    return spanVecs

def findFrame(spanVecs, bases):
    
    principle = []
    for vec in spanVecs + bases:
        
        if(rs.VectorLength(vec) <= EPSILON):
            continue
        
        included = isVectorPresent(principle, vec)
        
        if(not included):
            principle += [vec]
        
        if(len(principle) >= 2):
            # shortcut, x y vectors found
            break
    
    assert len(principle) == 2
    
    return principle

def isVectorPresent(pool, vec):
    
    included = False
    for poolVec in pool:
        angle = rs.VectorAngle(poolVec, vec)
        if(angle <= EPSILON):
            included = True
            break
            
    return included

def view(motion, isFreedom=False, bases=None):

    moType, axis, refPt, spanVecs = motion
    moType = moType if isFreedom else 'c'
    isTrans = moType == 't'
    
    # convert
    axis = [rs.AddPoint(vec) for vec in axis]
    refPt = rs.AddPoint(refPt)
    
    spanVecs = unpackSpans(spanVecs)
    
    if(bases == None):
        bases = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

    # create ordered list of reference frame vectors
    directors = axis + spanVecs
    xTar, yTar = findFrame(directors, bases)
    # create target frame using vectors
    targetFrame = rs.PlaneFromFrame(refPt, xTar, yTar)
    
    # find out which reference axial space to use
    if(not isTrans):
        axisDeg, spanDeg = len(axis), len(spanVecs)
    else:
        axisDeg, spanDeg = len(spanVecs), 0
    if(axisDeg == 1 and spanDeg == 1): # special case override
        axisDeg, spanDeg = 1, 0
    pathKey = treePath(axisDeg, spanDeg)
    axesToUse = refAxes.Branch(pathKey)
    spacesToUse = refSpaces.Branch(pathKey)
    frameToUse = refFrames.Branch(pathKey)[0]
    
    """
    if(axesToUse != None):
        vc.debug += axesToUse
    if(spacesToUse != None):
        vc.debug += spacesToUse
    vc.debug += [frameToUse] + [targetFrame]
    """

    # orient reference axial space to target
    linesOriented = []
    axesOriented = []
    spacesOriented = []
    for line in axesToUse:
        lineOriented = frameOrient(line, frameToUse, targetFrame, flags=1)
        linesOriented += [lineOriented]
        for obj in addAxesObj(lineOriented, moType):
            axesOriented += [obj]
    for obj in spacesToUse:
        objOriented = frameOrient(obj, frameToUse, targetFrame, flags=1)
        spacesOriented += [objOriented]
    
    return linesOriented, axesOriented, spacesOriented

def frameOrient(obj, source, target, flags=0):
    
    sourcePts = [source[0], source[0] + source[1], source[0] + source[2]]
    targetPts = [target[0], target[0] + target[1], target[0] + target[2]]
    
    oriented = rs.OrientObject(obj, sourcePts, targetPts, flags=flags)
    
    return oriented

def addAxesObj(line, moType):
    
    if(moType == 'r'):
        shapes = drawRotAxis(line, drawArc=True)
    elif(moType == 't'):
        shapes = drawTransAxis(line)
    elif(moType == 'c'):
        shapes = drawRotAxis(line, drawArc = False)
    
    return shapes

def drawArrow(anchor, vec, addCone=True):
    
    scale = rs.VectorLength(vec) / SCALE
    
    height = 2 * scale
    dia = 1 * scale
    
    end = rs.VectorAdd(anchor, vec)
    pipe = rs.AddPipe(rs.AddLine(anchor, end), 0, dia / 4, cap=1)
    if(addCone):
        apex = rs.VectorAdd(end, rs.VectorScale(rs.VectorUnitize(vec), height))
        plane = rs.PlaneFromNormal(apex, -vec)
        cone = rs.AddCone(plane, height, dia)
        arrow = rs.BooleanUnion(pipe + [cone])
    else:
        arrow = pipe
    
    return arrow[0]

def drawRotAxis(line, drawArc=True):
    
    scale = rs.CurveLength(line) / SCALE / 2
    
    height = 2 * scale
    dia = 1 * scale
    arcDeg = 240
    rotation = (360 - arcDeg) / 2
    arcRad = 2 * scale
    
    start = rs.CurveStartPoint(line) 
    mid = rs.CurveMidPoint(line)
    end = rs.CurveEndPoint(line)
    vecMS = rs.PointSubtract(start, mid)
    vecME = rs.PointSubtract(end, mid)
    vec = rs.VectorUnitize(rs.VectorSubtract(end, start))
    
    axisPipe = [drawArrow(mid, vecMS, False), drawArrow(mid, vecME, False)]

    #axisPipe = rs.AddPipe(line, 0, dia / 4, cap=2)
    arcPlane = rs.RotatePlane(rs.PlaneFromNormal(start, vec), rotation, vec)
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

def drawTransAxis(line):
    
    start = rs.CurveStartPoint(line) 
    mid = rs.CurveMidPoint(line)
    end = rs.CurveEndPoint(line)
    vecMS = rs.PointSubtract(start, mid)
    vecME = rs.PointSubtract(end, mid)
    
    arrows = [drawArrow(mid, vecMS), drawArrow(mid, vecME)]
    
    return arrows

def drawHighlightAxis(vec):
    
    if(dispCenter == None):
        center = rs.AddPoint(0, 0, 0)
    else:
        center = dispCenter
        
    ptStart = rs.PointAdd(center, rs.VectorScale(vec, SCALE))
    ptEnd = rs.PointAdd(center, rs.VectorScale(vec, -SCALE))
    line = rs.AddLine(ptStart, ptEnd)
    tube = rs.AddPipe(line, 0, .5, cap=1)

    return [line], tube

def drawTwins(cable, axes):
    
    cableLine = rs.AddLine(cable.start, cable.end)
    twins = [cableLine]
    for i in range(len(axes)):
        twins += rotateAboutAxis(twins, axes[i])
    
    original = twins.pop(0) # remove the original cable

    geo = []
    for line in twins:
        geo +=  rs.AddPipe(line, 0, cable.rad, cap=1)

    twins = [original] + twins

    return twins, geo

def rotateAboutAxis(objects, axis):

    if(dispCenter == None):
        center = rs.AddPoint(0, 0, 0)
    else:
        center = dispCenter
    
    result = rs.RotateObjects(objects, center, 180, axis, True)

    return result

def drawModes(info, bases):
    
    linesT = DataTree[object]()
    spacesT = DataTree[object]()
    visT = DataTree[object]()
    
    # draw mode-freedom spaces
    for i in range(len(info)):
        for j in range(len(info[i])):
            space = info[i][j]
            lines, axes, spaces = view(space, True, bases[i])
            
            path = treePath(i, j)
            linesT.AddRange(lines, path)
            visT.AddRange(axes, path)
            visT.AddRange(spaces, path)
            spacesT.AddRange(spaces, path)
    
    return linesT, spacesT, visT

def drawRods(info, bases):
    
    # rods and overconstraints

    rodsT = DataTree[object]()
    rodsLineT = DataTree[object]()
    rodsOCT = DataTree[object]()

    OCInfo = info[0]
    for rod in CM.model.rods:
        geo = rod.outputVis()

        line = rs.AddLine(rod.start, rod.end)
        rodsLineT.Add(line)
        
        if(rod.id in OCInfo):
            rodsOCT.AddRange(geo)
        else:
            rodsT.AddRange(geo)
    
    # spaces

    linesT = DataTree[object]()
    spacesT = DataTree[object]()
    visT = DataTree[object]()
    
    conSpaceInfo = info[1]
    for i in range(len(conSpaceInfo)):
        space = conSpaceInfo[i]
        lines, axes, spaces = view(space, False, bases)
            
        path = treePath(i)
        linesT.AddRange(lines, path)
        visT.AddRange(axes, path)
        visT.AddRange(spaces, path)
        spacesT.AddRange(spaces, path)

    return rodsLineT, rodsT, rodsOCT, linesT, spacesT, visT

def drawCables(info, bases):
    # cables
    OCInfo = info[0]
    twinInfo = info[3]

    cablesT = DataTree[object]()
    cablesOCT = DataTree[object]()
    cableTwinsT = DataTree[object]()
    twinsLineT = DataTree[object]()

    for i in range(len(CM.model.cables)):
        group = CM.model.cables[i]
        path = treePath(i)
        cablesT.AddRange([], path)
        cableTwinsT.AddRange([], path)
        twinsLineT.AddRange([], path)


        for cable in group:
            geo = cable.outputVis()
            if(cable.id in OCInfo[i]):
                cablesOCT.AddRange(geo)
            else:
                cablesT.AddRange(geo)
                axes = twinInfo[i]
                lines, geo = drawTwins(cable, axes)
                twinsLineT.AddRange(lines, path)
                cableTwinsT.AddRange(geo, path)



    # highlights
    hlInfo = info[1]
    highlightAxes = DataTree[object]()
    highlights = DataTree[object]()
    for i in range(len(hlInfo)):
        tempAxes = []
        tempGeo = []
        for j in range(len(hlInfo[i])):
            line, geo = drawHighlightAxis(hlInfo[i][j])
            tempAxes += line
            tempGeo += geo
        path = treePath(i)
        highlightAxes.AddRange(tempAxes, path)
        highlights.AddRange(tempGeo, path)

    # spaces
    conSpaceInfo = info[2]
    
    linesT = DataTree[object]()
    spacesT = DataTree[object]()
    visT = DataTree[object]()
    
    # draw group-constraint spaces
    for i in range(len(conSpaceInfo)):
        for j in range(len(conSpaceInfo[i])):
            space = conSpaceInfo[i][j]
            lines, axes, spaces = view(space, False, bases[i])
            
            path = treePath(i, j)
            linesT.AddRange(lines, path)
            visT.AddRange(axes, path)
            visT.AddRange(spaces, path)
            spacesT.AddRange(spaces, path)

    return cablesT, cablesOCT, cableTwinsT, twinsLineT, highlights, linesT, spacesT, visT

def drawSensors(info, bases):
    
    # sensors
    twinInfo = info[0]
    sensorsT = DataTree[object]()
    twinsLineT = DataTree[object]()
    sensorTwinsT = DataTree[object]()

    for i in range(len(CM.model.sensors)):
        group = CM.model.sensors[i]
        path = treePath(i)
        sensorsT.AddRange([], path)
        twinsLineT.AddRange([], path)
        sensorTwinsT.AddRange([], path)

        for sensor in group:
            geo = sensor.outputVis()
            sensorsT.AddRange(geo)
            axes = twinInfo[i]
            lines, geo = drawTwins(sensor, axes)
            twinsLineT.AddRange(lines, path)
            sensorTwinsT.AddRange(geo, path)

    # spaces
    conSpaceInfo = info[1]

    linesT = DataTree[object]()
    spacesT = DataTree[object]()
    visT = DataTree[object]()

    # draw group-constraint spaces
    for i in range(len(conSpaceInfo)):
        for j in range(len(conSpaceInfo[i])):
            space = conSpaceInfo[i][j]
            lines, axes, spaces = view(space, False, bases[i])
            
            path = treePath(i, j)
            linesT.AddRange(lines, path)
            visT.AddRange(axes, path)
            visT.AddRange(spaces, path)
            spacesT.AddRange(spaces, path)

    # icons
    iconInfo = info[2]
    
    iconsT = DataTree[object]()

    for i in range(len(iconInfo)):
        path = treePath(i)
        fileNames = [iconPath + name for name in iconInfo[i]]
        iconsT.AddRange(fileNames, path)

    return sensorsT, sensorTwinsT, twinsLineT, linesT, spacesT, visT, iconsT

if __name__ == "__main__":
    if(len(bases) != 0):
        bases = unpackBases(bases)
        
    if(len(info) != 0):
        _, _, step1Vis = drawModes(info[0], bases[0])
        rodsBake, rods, step2OC, step2DirGuide, step2PosGuide, step2ConSpace = drawRods(info[1], bases[1])
        cables, step3OC, cableTwins, cablesBake, step3Highlight, step3DirGuide, step3PosGuide, step3ConSpace = drawCables(info[2], bases[2])
        sensors, sensorTwins, sensorsBake, step4DirGuide, step4PosGuide, step4ConSpace, iconFiles = drawSensors(info[3], bases[3])