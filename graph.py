import numpy as np
import sympy as sp

class Graph(object):

    def __init__(self, dictForm):
        
        self.nodeCount = dictForm["nodeCount"]
        self.edgeCount = dictForm["edgeCount"]

        self.A = None
        self.C = None
        self.Q = None

        self._makeA(dictForm)
        self._makeC(dictForm)
        self._makeQ()
    
    def _makeA(self, dictForm):
        
        A = np.zeros((self.nodeCount, self.nodeCount)).astype(np.int)
        for id in dictForm["adjacency"]:
            e = dictForm["adjacency"][id]
            A[e[0], e[1]] = 1
            A[e[1], e[0]] = -1

        self.A = A
    
    def _makeC(self, dictForm):
        
        C = np.zeros((self.edgeCount, self.nodeCount))

        for id in dictForm["adjacency"]:
            e = dictForm["adjacency"][id]
            start, end = e[0], e[1]
            C[id, start] = -1
            C[id, end] = 1

        self.C = C
    
    def _makeQ(self):
        
        Q = sp.Matrix(self.C.T).nullspace()
        Q = np.asarray(Q).astype(np.float)
        
        self.Q = Q
    
    def _dijkstra(self, ground):

        nodes = set(range( self.nodeCount))

        dist = {}
        for i in range(self.nodeCount): dist[i] = float("inf")
        dist[ground] = 0
        prev = {}

        while(len(nodes) != 0):

            nodeDists = [(node, dist[node]) for node in nodes]
            cur = sorted(nodeDists, key=lambda x: x[1])[0][0]

            nodes.remove(cur)
            neighbors = np.nonzero(self.A[cur])[0]
            for n in neighbors:
                alt = dist[cur] + 1
                if(alt < dist[n]):
                    dist[n] = alt
                    prev[n] = cur

        return dist, prev

    def _makePathVec(self, path):

        pathVec = np.zeros(self.edgeCount)

        for i in range(len(path) - 1):
            nodeStart, nodeEnd = path[i], path[i + 1]
            en0, en1 = self.C[:, nodeStart], self.C[:, nodeEnd]
            e = np.nonzero(en0 * en1)[0][0]
            pathVec[e] = self.A[nodeStart, nodeEnd]
        
        return pathVec
    
    def makeP(self, ground, target):


        _, prev = self._dijkstra(ground)
        path = Graph._trackPath(prev, ground, target)
        
        P = self._makePathVec(path)

        return P
    
    @staticmethod    
    def _trackPath(prev, ground, target):

        cur = target
        path = [target]
        while(cur != ground):
            cur = prev[cur]
            path = [cur] + path
        
        return path
