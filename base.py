"""
This file provides basic functions to check for line segment intersection
and for checking the validity of a diagonal in a simple polygon
"""

from typing import List
from config import *
import sys
import os
import subprocess
import numpy as np


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
        n_ = [n[1], 0-n[0]]
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