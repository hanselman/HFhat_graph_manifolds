### This file demonstrates how to use the program tree_manifolds_HFhat.py with a
### few examples.

### To compute HFhat for a given weighted tree, enter the weights as a list, [weight at vertex 0, weight at vertex 1, ... ]
### and enter the edges as a list of pairs of indices
### Define a WeightedTree instance by WeightedTree(weights, edges)
### Call the function plumb on this weighted tree to get HFhat
### result.rank() gives the rank of HFhat
### result.homology() gives |H_1| (for rational homology spheres), which can be used to check if a manifold is an L-space


from tree_manifolds_HFhat import *

### compute HFhat for the manifold obtained by the E_8 plumbing (i.e. the Poincare homology sphere)
weights = [-2,-2,-2,-2,-2,-2,-2,-2]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5, 6), (4,7)]
E8 = WeightedTree(weights, edges)

result = plumb(E8)
print 'The rank of HFhat of the E8 plumbing is',
print result.rank()
print


### compute HFhat for the Brieskorn homology sphere \Sigma(3,5,7),
### using the plumbing tree in Figure 1 of Ozsvath-Szabo's "On the Floer homology of plumbed tree-manifolds"

weights = [-3,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5),
         (1,6), (6,7), (7,8), (8,9), (9,10), (10,11)]
brieskorn = WeightedTree(weights, edges)

result = plumb(brieskorn)
print 'For \Sigma(3,5,7),'
print 'the rank of HFhat is ',
print result.rank()
print 'the size of H_1 is ',
print result.homology()
print


### The following section of code computes rank of HFhat and the size of H_1 for a range weights
### on a given graph, and outputs the results in a file called output.txt
### It also makes note of the manifolds which are L-spaces
### The first column of output.txt lists the weights, the second column gives |H_1|, and the third gives rk(HFhat)
### if |H_1| = 0, then the manifold in question is not a rational homology sphere

edges = [(0, 1), (0,2), (0, 3)]
#### weights:
####
#### b
####   a  c   
#### d

f2 = open('output.txt', 'w')
f2.write('edges: [(0, 1), (0,2), (0, 3)] \n \n')

counter = 0
lspace_counter = 0
depth = 5
for a in range(-depth, 0):
    for b in range(-depth, 0):
        for c in range(-depth, 0):
            for d in range(-depth, 0):
                weights = [a,b,c,d]
                tree = WeightedTree(weights, edges)
                result = plumb(tree)

                counter += 1
                f2.write(str(weights) + '    ' + str(result.homology()) + '   ' + str(result.rank()))
                if result.rank() == result.homology():
                    f2.write('    ** L-space **')
                    lspace_counter += 1
                f2.write('\n')

f2.write('\n')
f2.write('number of examples: ' + str(counter))
f2.write('\n')
f2.write('number of L-spaces: ' + str(lspace_counter))
f2.close()
