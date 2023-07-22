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


def poly_vis_cgal_interior(verts: List[List[float]], query: [List[float]]) -> List[List[List[float]]]:
    """
    Computing visibility polygon of a query point (interior)
    :param verts: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise.
    :param query: Query point [x, y]
    :return: vertices of the visibility polygon
    """
    arg = poly_vis_path + ' '
    arg += '--verts='
    for v in verts:
        arg += str(v[0]) + ',' + str(v[1]) + ','
    arg = arg[:-1] + ' '
    arg += '--query=' + str(query[0]) + ',' + str(query[1]) + ' '
    arg += '--boundary=false'

    print('Running ' + arg)

    if os.name == 'posix':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True)
    if os.name == 'nt':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE)

    popen.wait()
    output = popen.stdout.read().decode("utf-8")
    output = output.split()
    output = [float(o) for o in output]

    vis_poly = []
    for i in range(0, len(output), 2):
        vis_poly.append([output[i], output[i+1]])

    return vis_poly


def poly_vis_cgal_boundary(verts: List[List[float]], query: [List[float]], pre_query: [List[float]]) -> List[List[List[float]]]:
    """
    Computing visibility polygon of a query point in a polygon with point on the boundary
    :param verts: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise.
    :param query: Query point [x, y]
    :param pre_query: Point on the boundary of the polygon preceding the query point
    :return: vertices of the visibility polygon
    """
    arg = poly_vis_path + ' '
    arg += '--verts='
    for v in verts:
        arg += str(v[0]) + ',' + str(v[1]) + ','
    arg = arg[:-1] + ' '
    arg += '--query=' + str(query[0]) + ',' + str(query[1]) + ' '
    arg += '--pre_query=' + str(pre_query[0]) + ',' + str(pre_query[1]) + ' '
    arg += '--boundary=true'

    if os.name == 'posix':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True)
    if os.name == 'nt':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE)

    popen.wait()
    output = popen.stdout.read().decode("utf-8")
    output = output.split()
    output = [float(o) for o in output]

    vis_poly = []
    for i in range(0, len(output), 2):
        vis_poly.append([output[i], output[i+1]])

    return vis_poly


def test_vis():
    #polygon = [[0, 0], [1, 1], [2, 0], [2, 4], [1, 2], [0, 4]]
    #query = [0.1, 0.3]
    #vis_polygon = poly_vis_cgal_interior(polygon, query)
    polygon = [[0, 0], [1, 1], [2, 0], [2, 4], [1, 2], [0, 4]]
    query = [1, 1]
    pre_query = [0, 0]
    vis_polygon = poly_vis_cgal_boundary(polygon, query, pre_query)


    vis_polygon = np.asarray(vis_polygon + [vis_polygon[0]])
    polygon = np.asarray(polygon + [polygon[0]])

    fig, ax = plt.subplots()
    plt.plot(polygon[:, 0], polygon[:, 1], 'k')
    plt.plot(vis_polygon[:, 0], vis_polygon[:, 1], 'b')
    plt.axis('equal')
    plt.show()
    plt.waitforbuttonpress()
    plt.close()


if __name__ == '__main__':
    test_vis()