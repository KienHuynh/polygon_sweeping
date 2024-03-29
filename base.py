"""
This file provides basic functions to check for line segment intersection
and for checking the validity of a diagonal in a simple polygon
"""

from typing import List
from config import *
import sys
import os
import struct
import subprocess
import numpy as np
from numpy.linalg import inv

from dijkstra import Graph, DijkstraSPF
from pslg import *


def area2(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])


def left(a, b, c):
    return area2(a, b, c) > 0


def left_on(a, b, c):
    return area2(a, b, c) >= 0


def colinear(a, b, c):
    if a == b or a == c or b == c:
        return True
    return area2(a, b, c) == 0


# def inCone(polygon, a, b):
#     n = len(polygon)
#     a0 = polygon[(polygon.index(a) - 1 + n) % n]
#     a1 = polygon[(polygon.index(a) + 1) % n]
#
#     if left_on(a, a1, a0):
#         return left(a, b, a0) and left(b, a, a1)
#
#     return not (left_on(a, b, a1) and left_on(b, a, a0))


def xor(a, b):
    return (a and not b) or (not a and b)


def between(a, b, c):
    if not colinear(a, b, c):
        return False;

    if a[0] != b[0]:
        return ((a[0] <= c[0]) and (c[0] <= b[0])) or ((a[0] >= c[0]) and (c[0] >= b[0]))
    else:
        return ((a[1] <= c[1]) and (c[1] <= b[1])) or ((a[1] >= c[1]) and (c[1] >= b[1]))


def intersectProp(a, b, c, d):
    if (colinear(a, b, c) or
            colinear(a, b, d) or
            colinear(c, d, a) or
            colinear(c, d, b)):
        return False

    return xor(left(a, b, c), left(a, b, d)) and xor(left(c, d, a), left(c, d, b))


def intersect(a, b, c, d):
    if intersectProp(a, b, c, d):
        return True
    elif between(a, b, c) or between(a, b, d) or between(c, d, a) or between(c, d, b):
        return True
    return False


def diagonalie(a, b, verts):
    c = verts[0]
    i = 0
    a = verts[a]
    b = verts[b]
    for i in range(len(verts)):
        c = verts[i]
        c1 = verts[(i+1) % len(verts)]
        if ((c != a) and (c1 != a) and (c != b) and (c1 != b) and intersect(a, b, c, c1)):
            return False

    return True


def in_cone(a, b, verts):
    a1 = verts[(a + 1) % len(verts)]
    a0 = verts[(a - 1) % len(verts)]
    a = verts[a]
    b = verts[b]

    if (left_on(a, a1, a0)):
        return left(a, b, a0) and left(b, a, a1)

    return not (left_on(a, b, a1) and left_on(b, a, a0))


def diagonal(a, b, verts):
    return in_cone(a, b, verts) and in_cone(b, a, verts) and diagonalie(a, b, verts)


def point_equal(x, y):
    """
    Check if two points are the same
    :param x: List of two numbers
    :param y: List of two numbers
    :return:
    """
    if x[0] == y[0] and x[1] == y[1]:
        return True
    return False


def create_vis_graph(vertices: List[List[float]], holes: List[List[List[float]]] = None) -> Graph:
    """
    Create a visibility graph based on the vertices of the polygon
    :param vertices: List[List[float]]
    :return: Graph object
    """

    graph = Graph()
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            vis = True
            for k in range(len(vertices)):
                k1 = k
                k2 = (k + 1) % (len(vertices))
                if (k1 == i) and (k2 == j):
                    vis = True
                    break

                if not diagonal(i, j, vertices) and ((i + 1) % len(vertices)) != j and ((j + 1) % len(vertices)) != i:
                    vis = False
                    break

                if intersectProp(vertices[i], vertices[j], vertices[k1], vertices[k2]):
                    vis = False
                    break

            if vis:
                d = l2_dist(vertices[i], vertices[j])
                graph.add_edge(i, j, d)
                graph.add_edge(j, i, d)
    return graph


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


def float2hex(num):
    return hex(struct.unpack('<Q', struct.pack('<d', num))[0])


def poly_vis_cgal_boundary(verts: List[List[float]], query: [List[float]], pre_query: [List[float]]) -> List[List[List[float]]]:
    """
    Computing visibility polygon of a query point in a polygon with point on the boundary
    :param verts: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise.
    :param query: Query point [x, y]
    :param pre_query: Point on the boundary of the polygon preceding the query point
    :return: vertices of the visibility polygon
    """
    arg = poly_vis_path + ' '
    verts_arg = ''
    for v in verts:
        v0_hex = float2hex(v[0])
        v1_hex = float2hex(v[1])
        verts_arg += '--verts ' + str(v0_hex) + ' --verts ' + str(v1_hex) + ' '
    arg += verts_arg
    q0_hex = float2hex(query[0])
    q1_hex = float2hex(query[1])
    arg += '--query ' + q0_hex + ' --query ' + q1_hex + ' '
    pre_q0_hex = float2hex(pre_query[0])
    pre_q1_hex = float2hex(pre_query[1])
    arg += '--pre_query ' + pre_q0_hex + ' --pre_query ' + pre_q1_hex + ' '
    arg += '--boundary=true'

    if os.name == 'posix':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True)
    if os.name == 'nt':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE)

    popen.wait()
    output = popen.stdout.read().decode("utf-8")
    output = output.split()
    try:
        output = [struct.unpack('!d', bytes.fromhex(str(o)))[0] for o in output]
    except ValueError as e:
        abc = 1

    vis_poly = []
    for i in range(0, len(output), 2):
        vis_poly.append([output[i], output[i+1]])

    return vis_poly


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


def find_seen_edge(p, d, polygon):
    nv = np.copy(d)
    nv[0], nv[1] = 0 - nv[1], nv[0]
    min_dist = 0
    edge = []
    edge_index = []

    for i in range(len(polygon)-1):
        p1 = polygon[i]
        p2 = polygon[i+1]
        if (above_below(p1, nv, p, True) * above_below(p2, nv, p, True)) < 0:
            nu = p1 - p2
            nu[0], nu[1] = 0 - nu[1], nu[0]
            u = p1
            nuv = np.stack((nv, nu))
            try:
                nuv = inv(nuv)
            except np.linalg.LinAlgError as E:
                continue
            else:
                pos_inv = nuv @ (np.asarray([p @ nv, u @ nu]).T)
                pos = pos_inv
                pos_deg = find_degen_intersection(p, nv, np.asarray([p1, p2]))
                if len(pos_deg) > 0:
                    if pos_deg @ nv - p @ nv == 0:
                        pos = pos_deg

            if len(edge) == 0:
                lowest = pos
                min_dist = get_norm(pos - p)
                edge = [p1, p2]
                edge_index = [i, i+1]
            else:
                dist = get_norm(pos - p)
                if dist < min_dist:
                    min_dist = dist
                    lowest = pos
                    edge = [p1, p2]
                    edge_index = [i, i+1]
    return edge, edge_index


def vis_decompose(polygon: List[List[float]], query: [List[float]], pre_query: [List[float]]) -> List[List[List[float]]]:
    """
    Compute the visibility region (VP) from the query point in the polygon, then decompose the original polygon
    into smaller subpolygons based on the visibility region: subtract VP from the polygon and we will get some smaller
    subpolygons remaining, the VP will be the parent polygon and the remaining subpolygons will be its children.
    The children are ordered counterclockwise around VP.
    :param verts: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise.
    :param query: Query point [x, y]
    :param pre_query: Point on the boundary of the polygon preceding the query point
    :return:
    """
    n = len(polygon)
    vis_polygon = poly_vis_cgal_boundary(polygon, query, pre_query)
    m = len(vis_polygon)
    v_head = [x for x in range(len(vis_polygon)) if point_equal(vis_polygon[x], query)]
    v_head = v_head[0]
    vis_polygon = vis_polygon[v_head:] + vis_polygon[:v_head]

    subpolygons = []
    i = 1 # Note that since the query point is a vertex of the polygon, it will always see the first edge and last edge
    j = 1
    while i != m - 1:
        if point_equal(vis_polygon[i + 1], polygon[j + 1]):
            i += 1
            j += 1
            continue

        # If there is a different with the next vertices, then there is a subpolygon
        subpoly = []

        # Special case where the visibility point is also a vertex of P
        ind = [x for x in range(len(polygon)) if point_equal(polygon[x], vis_polygon[i+1])]
        if len(ind) > 0:
            subpoly += polygon[j:(ind[0] + 1)]
            j = ind[0]
            i += 1
            subpoly = subpoly[-1:] + subpoly[:-1]
            subpolygons.append(subpoly)
            continue

        pi = np.asarray(vis_polygon[i + 1])
        direction = pi - vis_polygon[0]

        edge, edge_index = find_seen_edge(pi, direction, polygon[1:])
        edge_index[0] += 1
        edge_index[1] += 1

        if point_equal(edge[1], polygon[j+1]):
            subpoly.append(pi)
            r = [x for x in range(len(polygon)) if point_equal(polygon[x], vis_polygon[i+2])]
            r = r[0]
            subpoly += polygon[edge_index[1]:(r+1)]
            i += 2
            j = r
        else:
            subpoly += polygon[j:(edge_index[0] + 1)]
            subpoly.append(vis_polygon[i+1])
            i += 1
            j = edge_index[0]

        subpoly = subpoly[-1:] + subpoly[:-1]
        subpolygons.append(subpoly)

        #while not point_equal(vis_polygon[i], polygon[j + 1]):


    return (vis_polygon, subpolygons)


def hist_decompose(pslg: PSLG, a0):
    """
    Get the children of the histogram polygon
    :param pslg:
    :param a0:
    :return:
    """
    hist_polygon = [list(pslg.nodes[a0].pos)]
    subpolygons = []
    p = a0 + 1
    while p != a0:
        edge = []
        edges = [e for e in pslg.edges[p] if e.to.visited == False]

        if len(edges) == 0:
            break
        if len(edges) == 1:
            edge = edges[0]
        else:
            edge = [e for e in edges if e.type == 'vis']
            edge = edge[0] # This is the interface edge
            # Get the subpolygon to the right of this visibility edge
            subpoly = pslg.get_left_subpolygon(edge.to.id, edge.src.id)
            #subpoly = subpoly[1:] + subpoly[:1]
            subpolygons.append(subpoly)

        p = edge.to.id
        edge.src.visited = True
        hist_polygon.append(list(pslg.nodes[edge.src.id].pos))


    if (len(subpolygons) > 0):
        if point_equal(subpolygons[-1][0], pslg.nodes[a0].pos):
            subpolygons = subpolygons[:-1]
    if (len(subpolygons) > 0):
        if point_equal(subpolygons[0][1], pslg.nodes[a0+1].pos):
            subpolygons = subpolygons[1:]

    return (hist_polygon, subpolygons)


def above_below(p, n, a, use_normal=True):
    """
    Return a positive value if p is above the line xn - an = 0, negative otherwise
    If n is not a normal vector, it must be the tangent vector of the line
    :param p: numpy 2D vector, point for testing
    :param n: numpy 2D vector, normal vector
    :param a: numpy 2D vector
    :return:
    """
    if use_normal:
        return (p - a)@n #p@n - a@n
    else:
        n_ = [0-n[1], n[0]]
        return p@n_ - a@n_


def get_norm(v):
    """
    Calcualte the norm (i.e. length) of a vector, l2
    :param v: [x, y]
    :return: float value
    """
    return np.sqrt(v[0]**2 + v[1]**2)


def get_cos(u, v):
    """
    Get cosine distance between two vectors (in rad)
    :param u: [x, y]
    :param v: [x, y]
    :return:
    """
    product = u[0]*v[0] + u[1]*v[1]
    return product / (get_norm(u) * get_norm(v))


def is_reflex(a, b, c):
    """
    Check if the interior angle (a, b, c) is bigger than 180
    :param a: [x, y]
    :param b: [x, y]
    :param c: [x, y]
    :return:
    """
    ab = b - a
    ab[0], ab[1] = 0-ab[1], ab[0]
    above = above_below(c, ab, a)
    return above < 0


def is_vert_interior(a, b, c, d):
    """
    Check if the upward ray from b (upward = nv) is interior
    Only work if b is a reflex vertex :D
    :param a: [x, y]
    :param b: [x, y]
    :param c: [x, y]
    :param nv: [x, y]
    :return:
    """
    nv = np.copy(d)
    nv[0], nv[1] = 0 - nv[1], nv[0]
    above_a = above_below(a, nv, b, True)
    above_c = above_below(c, nv, b, True)

    if above_c == 0:
        if above_a > 0:
            return 1
        else:
            return -1

    if above_a == 0:
        # Upward ray align with ab
        if above_c < 0:
            return 1
        else:
            return -1


    if above_a * above_c > 0:
        return 1
    else:
        return -1


def l2_dist(p1, p2):
    """
    Compute l2 distance between p1 and p2
    :param p1: List[float]
    :param p2: List[float]
    :return:
    """
    p1 = np.asarray(p1).astype(float)
    p2 = np.asarray(p2).astype(float)
    return float(np.sqrt(np.sum((p1 - p2) ** 2)))