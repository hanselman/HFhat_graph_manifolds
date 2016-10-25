### This file demonstrates how to use the program graph_manifolds_HFhat.py with a
### few examples. This program is very similar to tree_manifolds_HFhat.py, with a
### few changes to allow for non-tree graphs with non-zero genus

### Graph manifolds are represented by a bi-weighted graph. The first sest of weights
### is the genus of a given vertex and the second set is the euler number

### To compute HFhat for a given graph manifold:
###     enter the genus as a list, [genus at vertex 0, genus at vertex 1, ... ]
###     enter the euler numbers as a list, [weight at vertex 0, weight at vertex 1, ... ]
###     enter the edges as a list of pairs of indices
###     enter edge signs as a list of integers (+1 or -1)
###     define a BiweightedTree instance by BiweightedTree(genus, weights, edges, edge_signs)
### Call the function plumb on this weighted tree to get HFhat
### result.rank() gives the rank of HFhat

### We can also compute the bordered invariants of a graph manifolds with boundary. To add an extra
### unfilled boundary at a given vertex of index i, add an edge of the form (i, '*')
### These false edges need to have an edge_sign assigned, though it does not matter what it is




from graph_manifolds_HFhat import *

### compute HFhat for the manifold obtained by the E_8 plumbing (i.e. the Poincare homology sphere)
genus = [0,0,0,0,0,0,0,0]
weights = [-2,-2,-2,-2,-2,-2,-2,-2]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5, 6), (4,7)]
edge_signs = [1,1,1,1,1,1,1]            #since this graph is a tree, the edge signs do not matter. We make them all +1 arbitrarily
E8 = BiweightedGraph(genus, weights, edges, edge_signs)

result = plumb(E8)
print 'The rank of HFhat of the E8 plumbing is',
print result.rank()
print



### the manifold T^3 = T^2 x S^1 can represented in several equivalent graphs
### we compute HF_hat using two of these
### (for their equivalence of these graphs, see Neumann, A Calculus for Plumbing, operation R5 in Proposition 2.1

print 'the rank of HFhat of T^3'
# Option 1: a single vertex with genus 1 and weight 0
genus = [1]
weights = [0]
edges = []
edge_signs = []
T3 = BiweightedGraph(genus, weights, edges, edge_signs)
result = plumb(T3)
print 'graph1:  ',
print result.rank()

# Option 2: two vertices with genus 0 and weight 0, connected by two edges of opposite sign
genus = [0,0]
weights = [0,0]
edges = [(0,1),(0,1)]
edge_signs = [1,-1]
T3 = BiweightedGraph(genus, weights, edges, edge_signs)
result = plumb(T3)
print 'graph2:  ',
print result.rank()
print


### similar to option 1 above, we can represent S_g x S^1, where S_g is the closed surface of genus g,
### by a graph with a single vertex and no edges. Here we compute HF_hat of these S^1 bundles for g up to 5
### note: if g<0, then S_g is a connect sum of g copes of RP^2. In this case, the bundle is not actually S_g x S^1,
###       but the trivial orientable bundle over S_g
print 'trivial orientable S^1 bundle over (surface of genus g)'
weights = [0]
edges = []
edge_signs = []
for g in range(-5,6):
    genus = [g]
    graph = BiweightedGraph(genus, weights, edges, edge_signs)
    print 'genus', g, ':  ',
    result = plumb(graph)
    print result.rank()
print


### Figure 6.5 in Gompf and Stipsicz, 4-manifolds and Kirby calculus, gives an example of a
### Kirby diagram for a 4-manifold obtained from plumbing, along with a corresponding plumbing graph
### the boundary of this 4-manifold is a graph manifold, represented by the same biweighted graph
### We can compute HF hat of this 3 manifold as follows:
genus = [1, 2, 0, -2]
weights = [3, -2, 4, 5]
edges = [(0,1),(0,2),(0,3)]
edge_signs = [-1,-1,-1]
T3 = BiweightedGraph(genus, weights, edges, edge_signs)
result = plumb(T3)
print 'Gompf-Stipsicz example (Figure 6.5):  ',
print result.rank()
print



### Graph manifolds with boundary:
### The following graph represents a manifold described by Liam Watson as a "Heegaard Floer solid tori"
### note that the edge (0,'*') indicates there is an unglued boundary component on the S^1 bundle
### corresponding to the vertex indexed by 0
print 'CFD of a graph manifold with one boundary component'
genus = [0,0,0,0,0]
weights = [0,2,-2,-2,2]
edges = [(0,1),(1,2),(0,3),(3,4),(0,'*')]
edge_signs = [1,1,1,1,1]
graph = BiweightedGraph(genus, weights, edges, edge_signs)
result = plumb(graph)
result.relabel_generators()
result.print_ops2()     #using print_ops2 instead of print_ops displays the operations sorted into connected components
