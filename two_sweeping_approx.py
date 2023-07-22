import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from dijkstra import Graph, DijkstraSPF
import copy

import subprocess
import sys
import os
import time

from typing import List
from base import *

poly_vis_path = './cgal_util/poly_vis_cgal/poly_vis_cgal'

def poly_vis_cgal(verts: List[List[float]], query: [List[float]]) -> List[List[List[float]]]:
    """
    Computing visibility polygon of a query point in a polygon
    :param verts: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise.
    :param query: Query point [x, y]
    :return: vertices of the visibility polygon
    """
    arg = poly_vis_path + ' '
    arg += '--verts='
    for v in verts:
        arg += str(v[0]) + ',' + str(v[1]) + ','
    arg = arg[:-1]
    arg += '--query=' + str(query[0]) + ',' + str(query[1])



if __name__ == '__main__':
    polygon = [[0, 0], [1, 1], [2, 0], [2, 4], [1, 2], [0, 4]]
    query = [0.1, 0.1]