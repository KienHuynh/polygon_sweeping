import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from numpy.linalg import inv
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



def find_degen_intersection(p, n, e):
    """
    Check if either endpoints of edge e belong to the line xn - pn
    :param p: [x, y]
    :param d: [x, y]
    :param e: [[x1, y1], [x2, y2]]
    :return:
    """
    d1 = e[0,:]@n - p@n
    d2 = e[1,:]@n - p@n
    if d1 == 0:
        if d2 == 0:
            return []
        else:
            return e[0, :]
    elif d2 == 0:
        return e[1, :]
    return []



def find_lowest_vert_intersection(p, d, pslg: PSLG):
    """
    Find the lowest intersection with pslg when shooting a ray in the d direction from p,
    return the intersection point and the edge in pslg the ray intersects with
    :param p: [x, y]
    :param d: [x, y]
    :param pslg: PSLG
    :return: ([x, y], Edge), [x, y] are the coordinates of the intersection and Edge is the edge in PSLG
    """
    nv = np.copy(d)
    nv[0], nv[1] = 0-nv[1], nv[0]
    lowest = []
    min_dist = 0
    edge = []
    for e in pslg.edge_set:
        if (e.src.pos[0] == p[0] and e.src.pos[1] == p[1]) or  (e.to.pos[0] == p[0] and e.to.pos[1] == p[1]):
            continue
        if e.type == "poly":
            #print(e.src.id, e.to.id)
            #print(e.src.pos, e.to.pos)
            if (above_below(e.src.pos, nv, p, True) * above_below(e.to.pos, nv, p, True)) <= 0:
                nu = e.src.pos - e.to.pos
                nu[0], nu[1] = 0-nu[1], nu[0]
                u = e.src.pos
                nuv = np.stack((nv, nu))
                try:
                    nuv = inv(nuv)
                except np.linalg.LinAlgError as E:
                    continue
                else:
                    pos_inv = nuv @ (np.asarray([p @ nv, u @ nu]).T)
                    pos = pos_inv
                    pos_deg = find_degen_intersection(p, nv, np.asarray([e.src.pos, e.to.pos]))
                    if len(pos_deg) > 0:
                        if pos_deg@nv - p@nv == 0:
                            pos = pos_deg

                if above_below(pos, d, p) <= 0:
                    continue

                if len(lowest) == 0:
                    lowest = pos
                    min_dist = get_norm(pos - p)
                    edge = e
                else:
                    dist = get_norm(pos - p)
                    if dist < min_dist:
                        min_dist = dist
                        lowest = pos
                        edge = e

    return lowest, edge


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
        right_polygon = np.vstack((right_polygon[1:, :], right_polygon[0, :].reshape(1,2)))
        right_polygon = np.asarray(poly_vis_cgal_boundary(right_polygon, right_polygon[0], right_polygon[-1]))

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

    fig, ax = plt.subplots()
    pslg.draw_adj_list()
    plt.plot([a[0], b[0]], [a[1], b[1]], 'r')
    plt.axis('equal')
    plt.show()
    plt.close()

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
        left_polygon = np.asarray(poly_vis_cgal_boundary(left_polygon, left_polygon[0], left_polygon[-1])) - 0.1

    fig, ax = plt.subplots()
    pslg.draw_adj_list()
    if len(right_polygon) > 0:
        draw_filled_polygon(right_polygon)
    if len(left_polygon) > 0:
       draw_filled_polygon(left_polygon)

    #plt.plot(left_polygon[[-1, 0], 0], left_polygon[[-1, 0], 1], 'g')
    plt.plot([a[0], b[0]], [a[1], b[1]], 'r')
    plt.axis('equal')
    plt.show()
    plt.close()

    # Retrieve the parent histogram polygon, as well as the children subpolygons
    ai = a0
    bi = b0

    #return main_hist