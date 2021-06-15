import argparse
import os.path, re, glob, sys, random
import matplotlib.pyplot as plt

"""
    tree_order_evaluation v1.0, 2021-06-15
    Reordering the internal nodes of a tree to obtain an ordering of the laves which
    minimise the number of conflicts with the lexicographic order of their labels.
    Copyright (C) 2020-2021 - Philippe Gambette, Olga Seminck
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
"""

# Parameter initialisation
currentFolder = os.path.abspath(os.path.dirname(sys.argv[0]))
inputFile = os.path.join(currentFolder,'input.txt')
outputFile = open(os.path.join(currentFolder,'output.txt'),"w")
testNb = 1000


"""
add all the leaves contained in clusterToPropagate to the value
corresponding to the key of all ancestors of node in dictionary cluster
"""
def propagateClusterAbove(node,clusterToPropagate,parent,cluster):
    #print("propagateClusterAbove("+str(node)+","+str(cluster)+")")
    if parent[node] > -1:
        if parent[node] not in cluster:
            cluster[parent[node]] = []
        for leaf in clusterToPropagate:
            cluster[parent[node]].append(leaf)        
        propagateClusterAbove(parent[node],clusterToPropagate,parent,cluster)


"""
parse line, a string at Newick format to build the tree structure
"""
def openNewick(line):
    i = 0
    currentDistanceFromRoot = 0
    currentNodeId = 0
    lastFoundNode = 0
    currentNode = ""
    currentCluster = []
    parent = {}
    label = {}
    parent[0] = -1
    cluster = {}
    mayBeALeaf = True
    hasASibling = False
    while i<len(line):
        if line[i] == "(":
            # update the distance from root for the new found node
            currentDistanceFromRoot += 1
            #print("Creating node "+str(lastFoundNode + 1)+"; parent: "+str(currentNodeId))
            previousNodeId = currentNodeId
            currentNodeId = lastFoundNode + 1
            lastFoundNode = currentNodeId
            # initialize an empty label for the current node
            currentNode = ""
            # initialize cluster (will contain the ids of the leaves found below the current node)
            currentCluster = []
            if not(currentNodeId in parent):
                parent[currentNodeId] = previousNodeId
                hasASibling = False                
            # print("Parent("+str(currentNodeId)+"): "+str(parent[currentNodeId]))
            mayBeALeaf = True            
        elif line[i] == ",":
            hasASibling = True
            # same distance from the root for the new found node
            # create an id for the new found node
            previousNodeId = currentNodeId
            currentNodeId = lastFoundNode + 1
            lastFoundNode = currentNodeId
            # show label of the previous node
            #print("Found node '"+currentNode+"'")
            # if previous node was a leaf: add it to the currentCluster to propagate:            
            if mayBeALeaf:
                res = re.search("^([^:]+):",currentNode)
                if res:
                    currentNode = res.group(1)
                label[previousNodeId] = currentNode            
                currentCluster.append(currentNodeId-1)
                # attribute the currently found cluster to the previous node
                cluster[currentNodeId-1] = currentCluster
                # propagate the currently found cluster to the ancestors:
                propagateClusterAbove(currentNodeId-1,currentCluster,parent,cluster)
            # initialize an empty label for the new found node
            currentNode = ""
            # reinitialize cluster
            currentCluster = []
            parent[currentNodeId] = parent[previousNodeId]
            #print("Parent2("+str(currentNodeId)+"): "+str(parent[currentNodeId]))
            mayBeALeaf = True
        elif line[i] == ")":
            # update the distance from root
            currentDistanceFromRoot -= 1
            # update the current node id: go back to the parent's id:
            previousNodeId = currentNodeId
            currentNodeId = parent[currentNodeId]
            # if previous node was a leaf: add it to the currentCluster to propagate:            
            if mayBeALeaf:
                res = re.search("^([^:]+):",currentNode)
                if res:
                    currentNode = res.group(1)
                label[previousNodeId] = currentNode            
                currentCluster.append(previousNodeId)
                # attribute the currently found cluster to the previous node
                cluster[previousNodeId] = currentCluster
                # propagate the currently found cluster to the ancestors:
                propagateClusterAbove(previousNodeId,currentCluster,parent,cluster)
            # show label of the previous node
            #print("Found node '"+currentNode+"'")
            currentNode = ""            
            mayBeALeaf = False
        else:
            currentNode += line[i]            
            
        i += 1
    
    children = {}
    for node in parent:
        if node > 0:
            if parent[node] not in children:
                children[parent[node]] = []
            children[parent[node]].append(node)
    
    
    result = {}
    result["parent"] = parent
    result["label"] = label
    result["cluster"] = cluster
    result["children"] = children
    
    #print("This tree was opened:")
    #print("Parents of vertices:"+str(parent))
    #print("Labels of vertices:"+str(label))
    #print("Clusters: "+str(cluster))
    #print ("Children: "+str(children))
    #print ("----------")
    return result


def allPermutationsAfterElement(list,i):
   result = []
   if i==len(list)-1:
      result = [list]
   else :
      for j in range(i,len(list)):
         list2 = list.copy()
         list2[i] = list[j]
         list2[j] = list[i]
         permutations = allPermutationsAfterElement(list2,i+1)
         for permutation in permutations:
            result.append(permutation)
   return result
      

def leafOrder(t,leaf):
   return t["label"][leaf]


def allPermutations(list):
   return allPermutationsAfterElement(list,0)


# count the number of crossings between the leaves of the clusters inside clusterList
# if all the clusters of clusterList are displayed in this order from left to right
def clusterCrossings(clusterList,t):
   crossingNb = 0
   # check all pairs of clusters
   for i in range(0,len(clusterList)):
      for j in range(i+1,len(clusterList)):
         clusterI = clusterList[i]
         clusterJ = clusterList[j]
         for leaf1 in clusterI:
            for leaf2 in clusterJ:
               if leafOrder(t,leaf1) > leafOrder(t,leaf2):
                  #print(str(leafOrder(t,leaf1))+">"+str(leafOrder(t,leaf2))+"!")
                  # we found a crossing!
                  crossingNb += 1
   #print("Cluster crossings in "+str(clusterList)+": "+str(crossingNb))
   return crossingNb


def permutationToStr(permutation,t):
   result = ""
   for vertex in permutation:
      result += "["
      for leaf in t["cluster"][vertex]:
         result += t["label"][leaf] + " "
      result += "]"
   return result

         
# recursively order the leaves of t,
# below node v, minimizing the number of crossings
# and return this number of crossings
def orderAndCountCrossingsBelow(t,v):
   #print("Visiting vertex "+str(v))
   crossings = 0
   if v not in t["children"]:
      #print(str(v)+" is a leaf")
      crossings = 0
   else:
      children = t["children"][v]
      #print(str(v)+" has "+str(len(children))+" children")
      permutations = allPermutations(children)
      #print("Try these orders of the children: "+str(permutations))
      bestPermutation = []
      bestCrossingNb = len(t["label"])
      for permutation in permutations:
         clusters = []
         for c in permutation:
            clusters.append(t["cluster"][c])
         crossingNb = clusterCrossings(clusters,t)
         #print("crossing number of vertex order "+str(permutation)+": "+str(crossingNb))
         if crossingNb < bestCrossingNb:
            bestCrossingNb = crossingNb
            bestPermutation = permutation.copy()
      crossings = bestCrossingNb
      t["children"][v] = bestPermutation
      #print("crossing number of best vertex order "+permutationToStr(bestPermutation,t)+": "+str(bestCrossingNb))
      for c in children:
         crossings += orderAndCountCrossingsBelow(t,c)
   return crossings
   

def printNewickTree(t):
   return printNewick(t,0)+";"


def printNewick(t,v):
   if v not in t["children"]:
      return str(t["label"][v])
   else:
      childNb = 0
      newick = ""
      for c in t["children"][v]:
         if childNb > 0:
            newick += ","
         newick += printNewick(t,c)
         childNb += 1
      return "("+newick+")"

# Create a pseudo-random order on n elements with Fisher-Yates shuffle
def pseudoRandomOrder(n):
   list = []
   for i in range(0,n):
      list.append(str(i))
   pseudoRandomList = []
   for i in range(0,n):
      pseudoRandomNumber = random.randint(i, n-1)
      pseudoRandomList.append(list[pseudoRandomNumber])
      list[pseudoRandomNumber] = list[i]
   return pseudoRandomList

with open(inputFile) as fd:
   i=0
   for line in fd.readlines():
      i += 1
      if i==1:
         # read the first line which should contain a tree in the Newick format
         t = openNewick(line)
   fd.close()
   # Reorder the tree while counting the number of crossings
   nbOfCrossings = orderAndCountCrossingsBelow(t,0)
   print("Correctly ordered dendrogram: "+printNewickTree(t))
   outputFile.writelines("Best ordered dendrogram:\n"+printNewickTree(t)+"\n")
   print("Nb of conflicts: "+str(nbOfCrossings))
   outputFile.writelines("Nb of conflicts: "+str(nbOfCrossings))
   
   outputFile.close()
   
   # Count the number of conflicts in the randomly re-ordered tree
   crossingNbList = []
   crossingNbDistribution = {}
   test = 0
   while test < testNb:
      order = pseudoRandomOrder(len(t["label"].keys()))
      test += 1
      #print("tree1"+str(printNewickTree(t)))
      newTree = openNewick(printNewickTree(t))
      i = 0
      
      for leaf in newTree["label"].keys():
         newTree["label"][leaf] = str(order[i])
         i += 1
      #print("tree2"+str(printNewickTree(newTree)))
      
      crossingNbNewTree = orderAndCountCrossingsBelow(newTree,0)
      crossingNbList.append(crossingNbNewTree)
      if crossingNbNewTree in crossingNbDistribution:
         crossingNbDistribution[crossingNbNewTree] += 1
      else:
         crossingNbDistribution[crossingNbNewTree] = 1
      
   minC = sorted(crossingNbDistribution.keys())[0]
   maxV = 0
   nbOfRandomOrderingsWithLessCrossings = 0
   for c in sorted(crossingNbDistribution.keys()):
      if c <= nbOfCrossings:
         nbOfRandomOrderingsWithLessCrossings += crossingNbDistribution[c]
      crossingNbDistribution[c] = crossingNbDistribution[c]*100/testNb
      if crossingNbDistribution[c] > maxV:
         maxV=crossingNbDistribution[c]
      print(str(c)+" conflicts: "+str(crossingNbDistribution[c])+"% of random orderings")
   print("For each nb of conflicts, percentage of random orderings having exactly this number of conflicts:")
   print(crossingNbDistribution)
   
   plt.bar(list(crossingNbDistribution.keys()), crossingNbDistribution.values(), color='g')
   plt.xlabel('Nb of conflicts between chronology and best dendrogram leaf ordering')
   plt.ylabel('Percentage of random orderings')
   plt.text(minC, maxV, r' '+str(nbOfCrossings)+' conflicts found for the best ordered dendrogram\n'+str(nbOfRandomOrderingsWithLessCrossings)+' random orderings over '+str(testNb)+' with at most '+str(nbOfCrossings)+' conflicts')
   plt.show()
