# Polygon Sweeping

This project deals with polygon sweeping using connected group(s) agents.

Requirement:
* Python 3.7.x
  * numpy
  * matplotlib
  * dijkstra

## ==========SWEEPING WITH ONLY TWO AGENTS==========

**Problem:** given a simple polygon, find a way to sweep it completely using two visibly connected agents. 
Each agent is represented as a point in 2D, both of them can be seen of a variable-length line segment.

Main file: two_sweeping.py

Requirement:
* CGAL 5.1: https://www.cgal.org/
* cxxopts (included in the project)

To run the python project:
1) Compile the files
```
./cgal/util/poly_decomp_cgal/poly_decomp_cgal.cpp
./cgal/util/poly_triangulation_cgal/poly_triangulation_cgal.cpp
```
using CGAL
2) Name them poly_decomp_cgal and poly_triangulation_cgal (with .exe for Windows)
3) Change the string variables that store the path to these executable files in two_sweeping.py
```
poly_decomp_path = './cgal_util/Debug/poly_decomp_cgal'
poly_triangulation_path = './cgal_util/Debug/poly_triangulation_cgal'
```
to wherever your actual executables are.

4) Run
```
python two_sweeping.py
```
with this line in the main function uncommented:
```
test_with_hole()
```
to test polygons with holes (or uncomment test_without_hole() in main for polygons without holes). Right now, you need to directly modify the polygons/holes in the functions test_with_hole and test_without_hole to test out different examples. The relevant variables are:

```
polygon = [[x1, y1], [x2, y2], ...]
holes = [
 [[x1, y1], [x2, y2], ...],
 [[x10, y10], [x11, y11], ...],
 ...
]
```
In which the [xi, yi] are **integer** coordinates of the vertices, written in counterclockwise order.
Then, pick any edge of the polygon and assign them as the starting location of the agents.
```
agents = [[x4, y4], [x5, y5]]
```
