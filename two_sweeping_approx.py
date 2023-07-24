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
    Return a positive value if p is below the line xn - an = 0, negative otherwise
    If n is not a normal vector, it must be the tangent vector of the line
    :param p: numpy 2D vector, point for testing
    :param n: numpy 2D vector, normal vector
    :param a: numpy 2D vector
    :return:
    """
    if use_normal:
        return p@n - a@n
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
        if e.type == "poly":
            if (above_below(e.src.pos, nv, p, True) * above_below(e.to.pos, nv, p, True)) < 0:
                nu = e.src.pos - e.to.pos
                nu[0], nu[1] = 0-nu[1], nu[0]
                u = e.src.pos
                nuv = np.stack((nv, nu))
                nuv = inv(nuv)
                pos = nuv@(np.asarray([p@nv, u@nu]).T)

                if len(lowest) == 0:
                    lowest = pos
                    min_dist = get_norm(pos - p)
                    edge = e
                else:
                    dist = get_norm(pos - p)
                    if dist < min_dist:
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

    # Check if b1 is below or above a0b0
    above = above_below(polygon[b1], nv, a, True)
    if above:
        cos_value = get_cos(b - a, polygon[b1] - polygon[a1])
        if cos_value >= 0:
            pos, edge = find_lowest_vert_intersection(b, nv, pslg)
            node = Node(len(pslg.nodes), pos)

            # Remove previous edge and add 4 new edges
            pslg.add_node(node)
            pslg.remove_edge(edge)
            pslg.add_edge(Edge(edge.src, node, 'poly'))
            pslg.add_edge(Edge(node, edge.to, 'poly'))
            pslg.add_edge(Edge(node, pslg.nodes[b0], 'vis'))
            pslg.add_edge(Edge(pslg.nodes[b0], node, 'vis'))

            fig, ax = plt.subplots()
            pslg.draw_adj_list()
            plt.axis('equal')
            plt.show()
            plt.waitforbuttonpress()
            plt.close()

            # fig, ax = plt.subplots()
            # draw_polygon(polygon)
            # plt.plot([b[0], pos[0]], [b[1], pos[1]], 'b')
            # plt.scatter([pos[0]], [pos[1]], color='r')
            # plt.axis('equal')
            # plt.show()
            # plt.waitforbuttonpress()
            # plt.close()


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
    edge1 = [[208, 560], [200, 552]]

    polygon2 = [[208, 560], [204, 564], [192, 560], [200, 564], [192, 572],
                [212, 572], [200, 584], [220, 580], [228, 584], [236, 576],
                [220, 568], [236, 556], [224, 560]]
    edge2 = [[208, 560], [204, 564]]

    polygon3 = [[208, 624],[224, 624],[236, 612],[232, 620],
            [252, 616],[248, 648],[216, 632],[232, 648],
            [212, 648],[224, 660],[188, 648],[212.249, 635.472],
            [192, 628],[199.586, 615.738],[208, 624]]
    edge3 = [[208, 624], [224, 624]]

    compute_histogram(polygon3, edge3)


if __name__ == '__main__':
    #test_vis()
    test_hist()