#!/usr/bin/env python
# encoding: utf-8
'''
 -- OPTICS - Density based clustering for spatial data

@author:     Aonghus Lawlor
            
@copyright:  2013. All rights reserved.
            
@license:   This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

            This program is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            GNU General Public License for more details.

            You should have received a copy of the GNU General Public License
            along with this program.  If not, see <http://www.gnu.org/licenses/>.

@contact:    aonghuslawlor@gmail.com

'''


import sys
import os
from rtree import index
from heapdict import heapdict
from priority_dict import priority_dict
from collections import namedtuple
import math
import functools
import itertools
import simplejson as json
import pickle
from optparse import OptionParser
import collections

import logging
logging.basicConfig(format="[%(asctime)s] {%(filename)s:%(funcName)s:%(lineno)d} %(levelname)s - %(message)s")
logging.getLogger().setLevel(logging.INFO)

LocalMax = collections.namedtuple('LocalMax', ['i', 'val'])

#@functools.total_ordering
class Node(object):
    def __init__(self, iStart, iEnd, parentNode):
        self.iStart = iStart
        self.iEnd = iEnd
        self.parentNode = parentNode
        self.range = self.iEnd - self.iStart
        self.iSplit = -1
        #self.children = set()
        self.children = []
        return

    def addChild(self, child):
        self.children.append(child)
        #self.children.add(child)
        return
    
    def split(self, i):
        self.iSplit = i
    
    def __str__(self):
        return "({}, {}, {}, {})".format(self.iStart, self.iEnd, self.iSplit, len(self.children))
    
    """
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.iStart == other.iStart and self.iEnd == other.iEnd
        else:
            return False

    def __lt__(self, other):
        if self.iStart == other.iStart:
            return self.iEnd < other.iEnd
        return self.iStart < other.iStart
    """
    
class OpticsClusterer(object):
    def __init__(self):
        return
    
    def checkLocalMax(self, i):
        j_low, j_high = max(0, i - self.localMaxWindow+1), min(self.nClusterOrder-1, i + self.localMaxWindow+1)
        reach_i = self.clusterOrder[i].reachability
        #logging.info("checkLocalMax: {} {} {} {}".format(i, j_low, j_high, reach_i))
        for j in xrange(j_low, j_high):
            if j != i and reach_i < self.clusterOrder[j].reachability:
                #logging.info("failed: {} {}".format(j, self.clusterOrder[j].reachability))
                return False
        return True
    
    def getLocalMaxima(self):
        localMaxima = []
        for i in xrange(1, self.nClusterOrder - 1):
            reach_i = self.clusterOrder[i].reachability
            
            if reach_i > self.clusterOrder[i-1].reachability and \
                reach_i >= self.clusterOrder[i+1].reachability and \
                self.checkLocalMax(i):                
                localMaxima.append(LocalMax(i=i, val=reach_i))
            #logging.info('{} {} {}'.format(i, reach_i, self.checkLocalMax(i)))
        return sorted(localMaxima, key=lambda x: x.val, reverse=False)
    
    def getAvReach(self, left=0, right=0):
            return float(sum(map(lambda x: self.clusterOrder[x].reachability, xrange(left, right)))) / (right - left)

    #@profile
    def clusterTree(self, node, parentNode=None, localmax=None):

        if localmax is None or len(localmax) == 0: #we're at a leaf...
            return 
        
        largestMax = localmax.pop()

        node.split(largestMax.i)

        #logging.info('clusterTree: node[{}] parentNode[{}] {} {} {}'.format(node, parentNode, largestMax.i, largestMax.val, len(localmax)))
        nodeL = Node(node.iStart, largestMax.i, node)
        nodeR = Node(largestMax.i + 1, node.iEnd, node)
        #localMaxL, localMaxR = [], []
        #for lm in localmax:
        #    if lm.i < largestMax.i: localMaxL.append(lm)
        #    if lm.i > largestMax.i: localMaxR.append(lm)
        
        nodeList = {'L': (nodeL, [lm for lm in localmax if lm.i < largestMax.i]), 
                    'R': (nodeR, [lm for lm in localmax if lm.i > largestMax.i]) }

        #logging.info('nodeL [{}] nodeR [{}]'.format(nodeList['L'][0], nodeList['R'][0]))
        
        if largestMax.val < self.reachabilityThreshold:
            node.split(-1)
            #logging.info('not significant- ignore and continue')
            self.clusterTree(node, parentNode, localmax)
            return

        checkL = int(self.checkFraction * nodeL.range)
        checkR = max(1, int(self.checkFraction * nodeR.range))
        avReachL = self.getAvReach(nodeL.iEnd - checkL, nodeL.iEnd)
        avReachR = self.getAvReach(nodeR.iStart, nodeR.iStart + checkR)

        reachLRatio = (avReachL / largestMax.val)
        reachRRatio = (avReachR / largestMax.val)

        #logging.info("ratio: {} {} {} : {} {}".format(largestMax.val, self.maxRatio, self.rejectionRatio, reachLRatio, reachRRatio))                 
        if reachLRatio > self.maxRatio or reachRRatio > self.maxRatio:
            #logging.info('testing...')
            if reachLRatio < self.rejectionRatio:
                #logging.info('removing R')
                del nodeList['R']
            if reachRRatio < self.rejectionRatio:
                #logging.info('removing L')
                del nodeList['L']
            if reachLRatio >= self.rejectionRatio and reachRRatio >= self.rejectionRatio:
                node.split(-1)
                #logging.info('ignore and continue {} {}'.format(node, parentNode))
                self.clusterTree(node, parentNode, localmax)
                return

        if nodeL.range < self.minClusterThreshold:
            if 'L' in nodeList: del nodeList['L']
        if nodeR.range < self.minClusterThreshold:
            if 'R' in nodeList: del nodeList['R']
        
        if len(nodeList) == 0:
            node.split(-1)
            return
        
        #logging.info('checking bypass')
        bypass = False
        if parentNode is not None:
            if float(node.range) / parentNode.range > self.similarityThreshold:
                #logging.info('parentNode: {} {} -> {}'.format(parentNode, node, parentNode.children))
                try:
                    parentNode.children.remove(node)
                except:
                    pass
                #parentNode.children.discard(node)
                bypass = True
        
        
        #for (k, (n, l)) in nodeList.iteritems():
        for k in ['L', 'R']:
            if k not in nodeList: continue
            (n, l) = nodeList[k]
            if bypass:
                parentNode.addChild(n)
                #logging.info('bypass {} {} {}'.format(k, n, parentNode))
                self.clusterTree(n, parentNode, l)
            else:
                node.addChild(n)
                #logging.info('no bypass {} {}'.format(n, node))
                self.clusterTree(n, node, l)
            
        return
    
    def getLeafNodes(self, node, leafNodes=[]):
        if node:
            if node.iSplit == -1:
                leafNodes.append(node)
            filter(lambda x : self.getLeafNodes(x, leafNodes), node.children)
        return False
    
    def run(self, clusterOrder=None):
        self.clusterOrder = clusterOrder
        self.nClusterOrder = len(self.clusterOrder)        

        #self.minPtsThreshold = 2
        self.similarityThreshold = 0.4
        self.localMaxFraction = 0.001                
        self.minClusterRatio = 0.005
        
        self.minClusterThreshold = min(5, int(self.minClusterRatio * self.nClusterOrder))
                        
        self.reachabilityThreshold = 0.003
        self.checkFraction = 0.8
        self.maxRatio = 0.75
        self.rejectionRatio = 0.7
        self.localMaxWindow = max(2, int(self.localMaxFraction * self.nClusterOrder))        
        
        #logging.info('self: {}'.format(str(self)))
        localmax = self.getLocalMaxima()        
        #logging.info('localmax {} {}'.format(len(localmax), localmax))
        #for l in localmax:
        #    logging.info('l: {}'.format(l))        
        
        self.root = Node(iStart=0, iEnd=self.nClusterOrder, parentNode=None)
        self.clusterTree(self.root, None, localmax)
        
        logging.info('root: {}'.format(self.root))
                        
        leafNodes = []
        self.getLeafNodes(self.root, leafNodes)
        
        clusters = [-1]  * self.nClusterOrder
        for c, leaf in enumerate(leafNodes):
            #print leaf
            for i in xrange(leaf.iStart, leaf.iEnd):
                #print X[clusterOrder[i].i]
                clusters[clusterOrder[i].i] = c
                #print clusters
        return clusters
    
    def __str__(self):
        return "nClusterOrder {}\nlocalMaxWindow {}\nminClusterThreshold {}\n".format(self.nClusterOrder, self.localMaxWindow, self.minClusterThreshold)
 
@functools.total_ordering
class ClusterEntry:
    def __init__(self, i=None, reachability=None):
        self.i = i
        self.reachability = reachability
        return
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.i == other.i
        else:
            return False

    def __lt__(self, other):
        if self.reachability == other.reachability:
            return self.i < other.i
        return self.reachability < other.reachability

    def __str__(self):
        #return "({}, {})".format(self.i, self.reachability)
        return str(self.__dict__)
    
    def __repr__(self):
        return self.__str__()
    
NNDist2 = namedtuple('NNDist2', ['id', 'distance'])

class Optics(object):
   
    def __init__(self, epsilon=sys.float_info.max, minPts=5):
        self.rtree = index.Rtree()
        self.processedIds = set()
        self.nPoints = None
        self.nodes = None
        self.epsilon = epsilon
        #self.epsilon2 = epsilon**2 if epsilon < sys.float_info.max else sys.float_info.max
        self.minPts = minPts
        #self.clusterOrder = Optics.ClusterOrder()
        self.clusterOrder = []
        return

    #@profile
    def getDistance(self, a, b): # for simple lists
        return math.sqrt( (a[0] - b[0]) ** 2 + (a[1] - b[1])**2 )
    
    #def getDistance2(self, a, b):
    #    return (a[0] - b[0]) ** 2 + (a[1] - b[1])**2
    
    #@profile
    def getNearestNeighboursDist(self, currentId):
        # get all nn's ordered by distance (minX, minY, maxX, maxY)
        pos = self.nodes[currentId]
        #return list(self.rtree.nearest( (pos[0], pos[1], pos[0], pos[1]), nPoints))
        res = []
        for nn in self.rtree.nearest( (pos[0], pos[1], pos[0], pos[1]), self.nPoints):
            #logging.info('nn: {} {}'.format(currentId, nn))
            if nn == currentId:
                continue
            # results are sorted by distance so we break when we get all with distance epsilon
            distance = self.getDistance(pos, self.nodes[nn])
            if distance <= self.epsilon:
                res.append( Optics.NNDist2(nn, distance) )
            else:
                break
        return res
    
    #@profile
    def getNearestNeighbours(self, currentId):
        # get all nn's ordered by distance (minX, minY, maxX, maxY)
        pos = self.nodes[currentId]
        return list(self.rtree.nearest( (pos[0], pos[1], pos[0], pos[1]), self.nPoints))


    #@profile  
    def expandClusterOrder(self, i):
        heap = heapdict()
        #heap = priority_dict()
        heap[i] = sys.float_info.max #float('inf')
        
        while heap:
            #logging.info('heap {}'.format( ','.join("({},{}) ".format(k,v) for k,v in heap.iteritems())))
            currentId, currentReachability = heap.popitem()            
            #currentId, currentReachability = heap.pop_smallest()
            
            #logging.info('processing: {} {} {}'.format(i, currentId, currentReachability))
            self.clusterOrder.append(ClusterEntry(i=currentId, reachability=currentReachability))
            if currentId in self.processedIds:
                logging.info('currentId in processedIds {}'.format(currentId))
            
            self.processedIds.add(currentId)
            
            neighbours = self.getNearestNeighbours(currentId)            
            n_neighbours = len(neighbours)

            if n_neighbours > 0 and n_neighbours > self.minPts:
                #coreDistance2 = neighbours[-1].distance2
                nn_last = neighbours[self.minPts]
                
                #coreDistance = neighbours[self.minPts-1].distance
                coreDistance = self.getDistance(self.nodes[currentId], self.nodes[nn_last])
                #logging.info('neighbours: {} {} -> {}'.format(n_neighbours, coreDistance2, neighbours))
                for nn in neighbours:
                    #logging.info('processedIds: {} {}'.format( (nn.id in self.processedIds), self.processedIds))
                    if nn == currentId or nn in self.processedIds: 
                        continue
                    
                    nn_distance = self.getDistance(self.nodes[currentId], self.nodes[nn])
                    newReachability = max(nn_distance, coreDistance)
                    oldReachability = heap.get(nn, None)
                    if oldReachability is None or newReachability < oldReachability:
                        heap[nn] = newReachability
        return

    #@profile
    def run(self, X):
        self.nodes = X
        self.nPoints = len(self.nodes)
        logging.info('#making rtree for {}...'.format(self.nPoints))
        for i, p in enumerate(self.nodes):
            self.rtree.add(i, tuple(p))
        logging.info('#made rtree: {}!'.format(self.rtree))
        
        for i, p in enumerate(self.nodes):
            if not i in self.processedIds:
                self.expandClusterOrder(i)
        return

    def opticsCluster(self, X):
        self.run(X)
        clusterer = OpticsClusterer()
        clusters = clusterer.run(self.clusterOrder)
        return clusters

def main():
    parser = OptionParser()
    parser.add_option("--infile", dest="infile", help="input file", default='optics_test.dat', type='str')
    parser.add_option("--npointscluster", dest="npointscluster", help="npointscluster", default=6, type='int')
    parser.add_option("--minPts", dest="minPts", help="minPts", default=6, type='int')
    parser.add_option("--epsilon", dest="epsilon", help="epsilon", default=sys.float_info.max, type='float')
    parser.add_option("--runtype", dest="runtype", help="runtype", default=0, type='int')

    (opts, args) = parser.parse_args()
    
    if opts.runtype == 0:
        X = []
        with open(opts.infile, 'r') as f:
            for line in f:
                if line.startswith('#') is False:
                    X.append(map(float, line.strip().split()[0:2]))
        #print 'X:', X
        optics = Optics(minPts=opts.minPts, epsilon=opts.epsilon)
        clusters = optics.opticsCluster(X)
        print clusters
    elif opts.runtype == 1:
        return
    elif opts.runtype == 2:
        return        
    return

PROFILE = 0
if __name__ == '__main__':
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'profile.txt'
        cProfile.run('main()', profile_filename)
        p = pstats.Stats(profile_filename)#, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative,sum')
        stats.print_stats()
        sys.exit(0)
    else:
        main()
