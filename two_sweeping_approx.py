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

poly_vis_path = './cgal_util/poly_vis_cgal/poly_vis_cgal'


def draw_polygon(polygon):
    n = len(polygon)
    plt.plot(polygon[:, 0], polygon[:, 1], 'k')
    plt.plot(polygon[[n-1, 0], 0], polygon[[n-1, 0], 1], 'k')


class Node:
    def __init__(self, _id, _pos):
        """
        2D position of the node
        :param _id:
        :param _pos:
        """
        self.id = _id
        self.pos = _pos
        self.edges = []
        self.visited = False


class Edge:
    def __init__(self, _from: Node, _to: Node, _type):
        """
        :param _from: Node object
        :param _to: Node onject
        :param _type: "poly" or "vis"
        """
        self.src = _from
        self.to = _to
        self.type = _type


class PSLG:
    def __init__(self, polygon):
        """

        :param polygon: Vertices of the polygon [[x1, y1], [x2, y2], ...], counter-clockwise
        """
        self.nodes = {}
        self.edges = {}
        self.edge_set = []
        n = len(polygon)
        for i in range(n):
            self.nodes[i] = Node(i, polygon[i])
        for i in range(n):
            edge = Edge(self.nodes[i], self.nodes[(i + 1) % n], "poly")
            self.edges[i] = [edge]
            self.edge_set.append(edge)

    def add_node(self, Node):
        num_node = len(self.nodes)
        self.nodes[num_node] = Node
        self.edges[num_node] = []


    def add_edge(self, _e):
        self.edge_set.append(_e)
        self.edges[_e.src.id].append(_e)


    def remove_edge(self, _e: Edge):
        """
        Remove edge _e from the edge set and the adjacency list
        :param _e:
        :return:
        """
        # Remove edge _e from the edge superset
        self.edge_set.remove(_e)

        # Remove edge _e from the adjacency list
        self.edges[_e.src.id].remove(_e)


    def draw_edge_set(self):
        for e in self.edge_set:
            c = 'k'
            if e.type == 'vis':
                c = 'b'

            src = e.src.pos
            to = e.to.pos
            plt.plot([src[0], to[0]], [src[1], to[1]], color=c)


    def draw_adj_list(self):
        for k, es in self.edges.items():
            for e in es:
                c = 'k'
                if e.type == 'vis':
                    c = 'b'

                src = e.src.pos
                to = e.to.pos
                plt.plot([src[0], to[0]], [src[1], to[1]], color=c)


    def node_exist(self, pos):
        for k, v in self.nodes:
            if v[0] == pos[0] and v[1] == pos[1]:
                return True

        return False


    def get_node(self, pos):
        for k, v in self.nodes.items():
            if v.pos[0] == pos[0] and v.pos[1] == pos[1]:
                return v

        return None


    def add_intersection(self, pos, edge, ray_start):
        # If this is an actual new node
        old_node = self.get_node(pos)
        if (old_node == None):
            node = Node(len(self.nodes), pos)

            # Remove previous edge and add 4 new edges
            self.add_node(node)
            self.remove_edge(edge)
            self.add_edge(Edge(edge.src, node, 'poly'))
            self.add_edge(Edge(node, edge.to, 'poly'))
            self.add_edge(Edge(node, self.nodes[ray_start], 'vis'))
            self.add_edge(Edge(self.nodes[ray_start], node, 'vis'))
            return node

        else:
            self.add_edge(Edge(old_node, self.nodes[ray_start], 'vis'))
            self.add_edge(Edge(self.nodes[ray_start], old_node, 'vis'))
            return old_node




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
                    pos = nuv @ (np.asarray([p @ nv, u @ nu]).T)

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
    above = above_below(polygon[b1], nv, a, True)
    if not (above > 0 and get_cos(b - a, polygon[b1] - polygon[a1]) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(b, nv, pslg)
        node = pslg.add_intersection(pos, edge, b0)

        # Ignore everything from b1 to node.id
        ai = node.id
        bi = edge.to.id

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
    if not (above >= 0 and get_cos(b - a, pslg.nodes[bi % n].pos - pslg.nodes[ai % n].pos) <= 0):
        # This means that the upward ray from b0 is interior to P and hits something
        # Find that intersection between the ray and the boundary
        pos, edge = find_lowest_vert_intersection(a, nv, pslg)
        node = pslg.add_intersection(pos, edge, a0)

    fig, ax = plt.subplots()
    pslg.draw_adj_list()
    plt.plot([a[0], b[0]], [a[1], b[1]], 'r')
    plt.axis('equal')
    plt.show()
    plt.close()

    # Retrieve the parent histogram polygon, as well as the children subpolygons



    #return main_hist

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


def test_hist():
    polygon1 = [[208, 560], [200, 552], [192, 560], [200, 564], [192, 572],
        [212, 572], [200, 584], [220, 580], [228, 584], [236, 576],
        [220, 568], [236, 556], [224, 560]]
    polygon1 = polygon1[::-1]
    edge1 = [[208, 560], [200, 552]]

    polygon2 = [[208, 560], [204, 564], [192, 560], [200, 564], [192, 572],
                [212, 572], [200, 584], [220, 580], [228, 584], [236, 576],
                [220, 568], [236, 556], [224, 560]]
    edge2 = [[208, 560], [204, 564]]

    polygon3 = [[208, 624],[224, 624],[236, 612],[232, 620],
            [252, 616],[248, 648],[216, 632],[232, 648],
            [212, 648],[224, 660],[188, 648],[212.249, 635.472],
            [192, 628],[199.586, 615.738]]
    edge3 = [[208, 624], [224, 624]]

    polygon4 = [[208, 624], [224, 624], [221.467, 629.45], [232, 620],
                [252, 616], [248, 648], [216, 632], [232, 648],
                [212, 648], [224, 660], [188, 648], [212.249, 635.472],
                [194.025, 628.747], [209.447, 628.807], [199.586, 615.738]]
    edge4 = [[199.586, 615.738], [208, 624]]

    polygon5 = [[240, 600], [240, 592], [248, 592], [248, 600], [252, 600],
        [252, 592], [256, 592], [256, 612], [252, 612], [252, 608],
        [248, 608], [248, 604], [244, 604], [244, 608], [236, 608], [236, 600]]

    polygon = polygon1
    for i in range(0, len(polygon)):
        edge_i = [polygon[i], polygon[(i + 1) % len(polygon)]]
        compute_histogram(polygon, edge_i)


if __name__ == '__main__':
    #test_vis()
    test_hist()