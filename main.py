from model import CompMech
from input import inputModel, viewCenter, ground, target, modeSwitch

motions = []
isFreedom = []
bases = []
result = []

if(modeSwitch != -1):
    model = CompMech(inputModel)
    DOFsys = model.findFreedom(ground, target)
    bases = str(DOFsys.outputBases())
    stageFrames = str(model.outputStageFrames())

if(modeSwitch == 0): # DOF analysis
    DOFsys.recenterSpaces(viewCenter)
    motions = str(DOFsys.output())
    isFreedom = True
elif(modeSwitch == 1): # constraint space
    conSys = DOFsys.constraintSys()
    conSys.recenterSpaces(viewCenter)
    motions = str(conSys.output())
    isFreedom = False
elif(modeSwitch == 2): # motion check
    result = DOFsys.isMotionPermissible(model.goals)
elif(modeSwitch == 3): # complementary motion check
    actuation = DOFsys.findSolution(model.goals[0])
    if(actuation is None): 
        result = None
    else:
        result = model.findMotion(ground, checkTarget, actuation)
elif(modeSwitch == 4): # show needed constrints
    isOverCon, overConsInfo = model.checkOverConstraint()
    motions, deltaInfo = DOFsys.constraintDelta(model.goals, viewCenter)
    motions = str(motions)
    result = ' '.join([overConsInfo, motions])
    isFreedom = False
elif(modeSwitch == 5): # show actuation space
    motions = DOFsys.viewInactSpace(viewCenter)
    motions = str(motions)
    isFreedom = False
elif(modeSwitch == 6): # check actuation
    result = model.checkActuations(ground, target)
    isFreedom = True
elif(modeSwitch == 7): #　single DOF-solve fast simulation
    result = model.tarSolveFastSim(ground, target)
elif(modeSwitch == 8): # multiple DOF fast simulation
    result = model.freeFastSim(ground, target)

print()
print(motions)
print(isFreedom)
print(bases)
print(result)