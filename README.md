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

To run the python project:
1) Compile the file ./poly_decomp_cgal/poly_decomp_cgal/main.cpp
2) Name it poly_decomp_cgal (or poly_decomp_cgal.exe for Windows)
3) Change the path to this executable file in two_sweeping.py
```
poly_decomp_path = './poly_decomp_cgal/Debug/poly_decomp_cgal'
```
to wherever your actual executable is.

4) Run
```
python two_sweeping.py
```
