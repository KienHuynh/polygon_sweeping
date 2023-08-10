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
from animate import *



def compute_hist_vist_polygon(polygon, e):
    """
    Compute the histogram polygon from base edge e
    No recursion here
    :param polygon:
    :param polygon: Vertices of the polygon [[x1, y1], [x2, y2], ...], counterclockwise.
    :param e: Edge [[x1, y1], [x2, y2]], must be an edge of the polygon, counterclockwise as well.
    :return:
    """
    polygon = [list(p) for p in polygon]
    e = [list(i) for i in e]
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
    child_subpolygons = []

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
    right_vis_polygon = []
    right_children = []
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
        # right_polygon = np.vstack((right_polygon[1:, :], right_polygon[0, :].reshape(1,2)))
        right_polygon = right_polygon[1:] + right_polygon[:1]
        right_vis_polygon, right_children = vis_decompose(right_polygon, right_polygon[0], right_polygon[-1])

    # Start going around P, compute the projected image of each reflex vertex onto the boundary of P with direction "up"
    while bi != a0:
        bj = (bi + 1) % n
        if (is_reflex(pslg.nodes[ai].pos, pslg.nodes[bi].pos, pslg.nodes[bj].pos)):
            interior_chord = is_vert_interior(pslg.nodes[ai].pos, pslg.nodes[bi].pos, pslg.nodes[bj].pos, nv)
            if interior_chord == 1:
                try:
                    pos, edge = find_lowest_vert_intersection(pslg.nodes[bi].pos, nv, pslg)
                    node = pslg.add_intersection(pos, edge, bi)
                except IndexError as x:

                    abc = 1

                # if (edge.to.id - a0) % n < (bi - a0) % n:
                #     ai = bi
                #     bi = bj
                # else:
                #     ai = node.id
                #     bi = edge.to.id
                # continue

        ai = bi
        bi = bj

    # Take special care of the last edge
    bj = b0
    # Check if ai is below or above a0b0
    above = above_below(pslg.nodes[ai].pos, nv, a, True)
    left_vis_polygon = []
    left_children = []
    if not (above > 0 and get_cos(b - a, pslg.nodes[bi].pos - pslg.nodes[ai].pos) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(a, nv, pslg)
        node = pslg.add_intersection(pos, edge, a0)

        left_polygon = pslg.get_left_subpolygon(a0, node.id)
        left_vis_polygon, left_children = vis_decompose(left_polygon, left_polygon[0], left_polygon[-1])

    hist_polygon, hist_children = hist_decompose(pslg, a0)

    return left_vis_polygon, left_children, hist_polygon, hist_children, right_vis_polygon, right_children, right_polygon



def compute_hist_vis_polygoon_re(polygon, e, r2l):
    """
    Compute histogram + left/right vis polygon(s) from base edge e,
    then recursively compute that of casted shadow edges
    :param polygon: Vertices of the polygon [[x1, y1], [x2, y2], ...], counterclockwise.
    :param e: Edge [[x1, y1], [x2, y2]], must be an edge of the polygon, counterclockwise as well.
    :return:
    """
    polygon = [list(p) for p in polygon]
    e = [list(i) for i in e]
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
    child_subpolygons = []

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
    right_vis_polygon = []
    right_children = []
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
        right_vis_polygon, right_children = vis_decompose(right_polygon, right_polygon[0], right_polygon[-1])

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

                # if (edge.to.id - a0) % n < (bi - a0) % n:
                #     ai = bi
                #     bi = bj
                # else:
                #     ai = node.id
                #     bi = edge.to.id
                # continue

        ai = bi
        bi = bj

    # Take special care of the last edge
    bj = b0
    # Check if ai is below or above a0b0
    above = above_below(pslg.nodes[ai].pos, nv, a, True)
    left_polygon = []
    left_vis_polygon = []
    left_children = []
    if not (above > 0 and get_cos(b - a, pslg.nodes[bi].pos - pslg.nodes[ai].pos) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(a, nv, pslg)
        node = pslg.add_intersection(pos, edge, a0)

        left_polygon = pslg.get_left_subpolygon(a0, node.id)
        left_vis_polygon, left_children = vis_decompose(left_polygon, left_polygon[0], left_polygon[-1])


    hist_polygon, hist_children = hist_decompose(pslg, a0)

    child_subpolygons += left_children[::-1] + hist_children[::-1] + right_children[::-1]

    # Retrieve the parent histogram polygon, as well as the children subpolygons
    ai = a0
    bi = b0

    # Merge the left/right vis polygons and the histogram polygon into one
    main_node = HistogramNode(r2l=r2l)
    main_node.dividers = []
    main_node.subpolygons = []
    if (len(right_vis_polygon) > 0):
        dividing_edge = [right_vis_polygon[0], right_vis_polygon[-1]]
        left_edge = [right_vis_polygon[0], right_vis_polygon[-1]]
        right_edge = [right_vis_polygon[0], right_vis_polygon[1]]
        subpolygon = SubPolygon(right_vis_polygon, 'vis', None, left_edge, right_edge)
        main_node.dividers.append(dividing_edge)
        main_node.subpolygons.append(subpolygon)

    base_edge = [hist_polygon[0], hist_polygon[1]]
    left_edge = [hist_polygon[0], hist_polygon[-1]]
    right_edge = [hist_polygon[1], hist_polygon[2]]
    subpolygon = SubPolygon(hist_polygon, 'hist', base_edge, left_edge, right_edge)
    main_node.subpolygons.append(subpolygon)

    if (len(left_vis_polygon) > 0):
        dividing_edge = [left_vis_polygon[0], left_vis_polygon[1]]
        left_edge = [left_vis_polygon[0], left_vis_polygon[-1]]
        right_edge = [left_vis_polygon[0], left_vis_polygon[1]]
        subpolygon = SubPolygon(left_vis_polygon, 'vis', None, left_edge, right_edge)
        main_node.dividers.append(dividing_edge)
        main_node.subpolygons.append(subpolygon)

    main_node.compute_outer_boundary()

    # if r2l:
    #     draw_filled_polygon(main_node.outer_boundary, 'r')
    # else:
    #     draw_filled_polygon(main_node.outer_boundary, 'b')

    r2l = not r2l
    for c in child_subpolygons:

        c_node = compute_hist_vis_polygoon_re(c, [c[0], c[1]], r2l=r2l)
        main_node.children.append(c_node)

    main_node.interface_edge = [polygon[a0], polygon[b0]]

    #main_node.interface_edge

    return main_node


def get_max_reflex_chain(polygon, i):
    """
    Given a vertex of the polygon, find the longest reflex chain starting from i going counter-clockwise
    :param polygon:
    :param i:
    :return:
    """
    n = len(polygon)
    chain = [polygon[i]]
    #length = l2_dist(chain[0], chain[1])
    length = 0
    polygon = np.asarray(polygon)
    a0 = polygon[i % n, :]
    b0 = polygon[(i + 1) % n, :]
    vec0 = b0 - a0

    for j in range(i, n + i):
        ja = j % n
        jb = (j + 1) % n
        a1 = polygon[ja, :]
        b1 = polygon[jb, :]
        if j == i or above_below(b1, vec0, a0, use_normal=False) <= 0:
            chain.append(list(b1))
            vec0 = b1 - a1
            a0 = a1
            length += l2_dist(b1, a1)
        else:
            break

    return chain, length


def longest_reflex_chain(polygon: List[List[float]]):
    """
    Find the longest reflex chain in a polygon
    :param polygon: [[x1, y1], [x2, y2], ...]
    :return:
    """
    n = len(polygon)
    max_length = 0
    max_chain = []
    for i in range(n):
        chain, length = get_max_reflex_chain(polygon, i)
        if length > max_length:
            max_length = length
            max_chain = chain

    return max_chain


def lrc_decompose(polygon):
    chain = longest_reflex_chain(polygon)
    m = len(chain)

    main_node = HistogramNode(r2l=True)
    main_node.dividers = []
    main_node.subpolygons = []
    child_subpolygons = []

    next_polygon = polygon
    # Recall that the "main node" is made of smaller alternative subpolygons of type "vis" and "hist"
    # If we consider the ordering of the reflex chain (which is counter clockwise) to be left-to-right
    # these smaller subpolygons will be ordered from "right to left", i.e. when sweeping these subpolygons,
    # one endpoint will move clockwise along the reflex chain, but counter clockwise along the top chain
    # BUT the children will be ordered clockwise

    for i in range(m-1):
        a = chain[i]
        b = chain[(i+1) % m]
        e = [a, b]
        left_vis_polygon, left_children, hist_polygon, hist_children, right_vis_polygon, right_children, next_polygon \
            = compute_hist_vist_polygon(next_polygon, e)

        if (len(left_vis_polygon) > 0):
            dividing_edge = [left_vis_polygon[0], left_vis_polygon[1]]
            left_edge = [left_vis_polygon[0], left_vis_polygon[-1]]
            right_edge = [left_vis_polygon[0], left_vis_polygon[1]]
            subpolygon = SubPolygon(left_vis_polygon, 'vis', None, left_edge, right_edge)
            main_node.dividers.append(dividing_edge)
            main_node.subpolygons.append(subpolygon)

            child_subpolygons += left_children[::-1]

        base_edge = [hist_polygon[0], hist_polygon[1]]
        left_edge = [hist_polygon[0], hist_polygon[-1]]
        right_edge = [hist_polygon[1], hist_polygon[2]]
        subpolygon = SubPolygon(hist_polygon, 'hist', base_edge, left_edge, right_edge)
        main_node.subpolygons.append(subpolygon)

        child_subpolygons += hist_children[::-1]

        if i == m-2 and len(right_vis_polygon) > 0:
            dividing_edge = [right_vis_polygon[0], right_vis_polygon[-1]]
            left_edge = [right_vis_polygon[0], right_vis_polygon[-1]]
            right_edge = [right_vis_polygon[0], right_vis_polygon[1]]
            subpolygon = SubPolygon(right_vis_polygon, 'vis', None, left_edge, right_edge)
            main_node.dividers.append(dividing_edge)
            main_node.subpolygons.append(subpolygon)

            child_subpolygons += right_children[::-1]

    main_node.subpolygons = main_node.subpolygons[::-1]
    main_node.dividers = main_node.dividers[::-1]

    main_node.compute_outer_boundary()
    for c in child_subpolygons:
        c_node = compute_hist_vis_polygoon_re(c, [c[0], c[1]], r2l=False)
        main_node.children.append(c_node)

    return main_node

