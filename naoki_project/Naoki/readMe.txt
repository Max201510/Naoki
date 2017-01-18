PointCor.java: to describe a point coordinate
 value: the coordinate in one-dimensional space
 x: x value in two-dimensional space
 y: y value in two-dimensional space
 
Edge.java: to describe an edge information.
 index: edge index 
 direction: 1 means horizontal, 0 means vertical, 2 means inclined, 3 means arbitary
 length: edge length
 cindex: the connected component index to which this edge belongs

CC.java: a connected component -  a collection of edges
 cindex: connected componet index
 edges: an ArrayList of edges
 length: this connected componet lenght - add all edges lengths

Boundary.java: connected component representation - a collection of double numbers
 bound: an ArrayList of double numbers

Naoki.java: 
 boundTheshold: it is a parameter that decides to filter out some small boundaries
 folders: the connected componet folder
