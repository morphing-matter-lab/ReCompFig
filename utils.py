import numpy as np
import sympy as sp
from parameters import EPSILON, PREC, EUCSPACEDIM
import math as m
import decimal

decimal.getcontext().prec = PREC
spatialFrame = np.concatenate([np.zeros((1, 3)), np.identity(3)])

def _fixFloat(inp):
    if(isinstance(inp, np.ndarray)):

        inp[np.isclose(inp, np.zeros(inp.shape), atol=EPSILON)] = 0
        inp[np.isclose(inp, np.ones(inp.shape), atol=EPSILON)] = 1
    else:
        if(abs(inp) <= EPSILON): inp = 0
        elif(abs(1 - abs(inp)) <= EPSILON): inp = -1 if inp < 0 else 1

    diff = inp - _fixFloat(inp) 
    if(isinstance(diff, np.ndarray)):
        if(np.any(diff > 0)):
            print(diff)
    elif(diff > 0):
        print(diff)
    return inp

def fixFloat(inp):
    # recursively fix an array
    if(isinstance(inp, np.ndarray) or\
       isinstance(inp, list) or isinstance(inp, tuple)):
        converted = []
        for item in inp:
            converted += [fixFloat(item)]
        
        # 
        if(isinstance(inp, np.ndarray)):
            return np.asarray(converted)
        elif(isinstance(inp, tuple)):
            return tuple(converted)
        else:
            return converted
    else:
        fixed = round(inp, PREC + 1)
        return fixed

def p2t(point):
    # convert a point coordinate into a tuple
        
    return tuple([num for num in point])

def solve(eqs, symbols):

    # try to solve symbolically
    solution = sp.solve(eqs)
    solved = not isinstance(solution, list)
    
    if(solved):
        symbVals = [solution[symb] for symb in symbols]
        symbVals = np.asarray(symbVals).astype(np.float)
        return symbVals
    
    # try to solve numerically
    initVal = np.random.rand(len(symbols))
    for i in [9, 8, 7, 6, 5, 4]:
        try:
            solution = sp.nsolve(eqs, symbols, initVal, prec=i)
            solution = np.asarray(solution).reshape(-1).astype(np.float)
        except:
            continue
    
    assert len(solution) == len(symbols), "multiple solutions found"
    symbVals = solution
    
    return symbVals

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

def exhaustiveSubclasses(targClass):
    result = [targClass]

    subclasses = targClass.__subclasses__()
    if(len(subclasses) != 0):
        for subclass in subclasses:
            moreSubclassses = exhaustiveSubclasses(subclass)
            result += moreSubclassses

    return result
