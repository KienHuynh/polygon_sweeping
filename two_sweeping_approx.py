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
from pslg import *
from config import *
from graph import *
from animate import *


def compute_histogram(polygon, e):
    """
    Compute histogram polygon from base edge e
    :param polygon: Vertices of the polygon [[x1, y1], [x2, y2], ...], counterclockwise.
    :param e: Edge [[x1, y1], [x2, y2]], must be an edge of the polygon, counterclockwise as well.
    :return:
    """
    a = e[0]
    b = e[1]
    a0 = polygon.index(a)
    b0 = polygon.index(b)
    n = len(polygon)
    assert (b0 - a0) % n == 1, "Given e is not an edge of polygon"

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    polygon = np.asarray(polygon, dtype=np.float64)

    # Create a PSLG out of the polygon
    pslg = PSLG(polygon)

    # Get normal vector of edge a0b0, rotated counter clockwise
    nv = b - a
    nv[0], nv[1] = 0 - nv[1], nv[0]

    # Take special care of the first edge
    a1 = b0
    b1 = (b0 + 1) % n
    ai = b0
    bi = b1

    # Check if b1 is below or above a0b0
    right_polygon = []
    above = above_below(polygon[b1], nv, a, True)
    if not (above > 0 and get_cos(b - a, polygon[b1] - polygon[a1]) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(b, nv, pslg)
        node = pslg.add_intersection(pos, edge, b0)

        # Ignore everything from b1 to node.id
        ai = node.id
        bi = edge.to.id

        # Compute the visibility polygon from b0, restricted to the right of b0 and node
        right_polygon = pslg.get_left_subpolygon(ai, b0)
        # Shift the polygon by 1 because right now the first vertex is not b0
        #right_polygon = np.vstack((right_polygon[1:, :], right_polygon[0, :].reshape(1,2)))
        right_polygon = right_polygon[1:] + right_polygon[:1]
        right_vis_polygon, right_subpolygons = vis_decompose(right_polygon, right_polygon[0], right_polygon[-1])

        # Decompose right_polygon into smaller subpolygons: the vis polygon will be the parent, subtract the main polygon
        # with the right vis polygon and we will have its children


    # Start going around P, compute the projected image of each reflex vertex onto the boundary of P with direction "up"
    while bi != a0:
        bj = (bi + 1) % n
        if (is_reflex(pslg.nodes[ai].pos, pslg.nodes[bi].pos, pslg.nodes[bj].pos)):
            interior_chord = is_vert_interior(pslg.nodes[ai].pos, pslg.nodes[bi].pos, pslg.nodes[bj].pos, nv)
            if interior_chord == 1:
                pos, edge = find_lowest_vert_intersection(pslg.nodes[bi].pos, nv, pslg)
                node = pslg.add_intersection(pos, edge, bi)

                if (edge.to.id - a0) % n < (bi - a0) % n:
                    ai = bi
                    bi = bj
                else:
                    ai = node.id
                    bi = edge.to.id
                continue

        ai = bi
        bi = bj

    # Take special care of the last edge
    bj = b0
    # Check if ai is below or above a0b0
    above = above_below(pslg.nodes[ai].pos, nv, a, True)
    left_polygon = []
    if not (above > 0 and get_cos(b - a, pslg.nodes[bi].pos - pslg.nodes[ai].pos) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(a, nv, pslg)
        node = pslg.add_intersection(pos, edge, a0)

        left_polygon = pslg.get_left_subpolygon(a0, node.id)
        left_vis_polygon, left_subpolygons = vis_decompose(left_polygon, left_polygon[0], left_polygon[-1])

    fig, ax = plt.subplots()
    pslg.draw_adj_list()
    if len(right_polygon) > 0:
        draw_filled_polygon(right_vis_polygon)
        for p in right_subpolygons:
            draw_polygon(p, 'r')

    if len(left_polygon) > 0:
       draw_filled_polygon(left_vis_polygon)
       for p in left_subpolygons:
           draw_polygon(p, 'r')
       draw_polygon(left_polygon, 'g')

    #plt.plot(left_polygon[[-1, 0], 0], left_polygon[[-1, 0], 1], 'g')
    plt.plot([a[0], b[0]], [a[1], b[1]], 'r')
    plt.axis('equal')
    plt.show()
    plt.close()

    # Retrieve the parent histogram polygon, as well as the children subpolygons
    ai = a0
    bi = b0

    #return main_hist