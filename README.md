pbfs
====

Parallel Breadth First Search

The code implement BFS in parallel. It uses Intel Cilk++ library and runs on share memory system. 
A data structure named “Bag of Bone” is applied instead of Queue.This data structure can be refered to this paper: http://dl.acm.org/citation.cfm?id=1810534

The program can traverse a graph with more than 2^25 vertices and edges within a second.
