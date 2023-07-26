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
from config import *

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


def is_same_poly(polys, edge):
    for poly in polys:
        if (edge[0] in poly) and (edge[1] in poly):
            return True

    return False


def is_actual_edge(polys, edge):
    for poly in polys:
        n = len(poly)
        if (edge[0] in poly) and (edge[1] in poly):
            i = poly.index(edge[0])
            j = poly.index(edge[1])
            if ((i + 1) % n) == j or ((j + 1) % n) == i:
                return True

    return False


def create_vis_graph_with_holes(vertices: List[List[float]], holes: List[List[List[float]]]) -> Graph:
    """
    Create a visibility graph based on the vertices of the polygon with holes (can be empty)
    :param vertices: List[List[float]]
    :return: Graph object
    """

    graph = Graph()
    vertices = copy.deepcopy(vertices)
    main_poly = copy.deepcopy(vertices)
    all_polys = [main_poly] + holes

    for h in holes:
        vertices += h

    for h in holes + [main_poly]:
        for i in range(len(h)):
            j = (i + 1) % len(h)
            d = l2_dist(h[i], h[j])
            ii = vertices.index(h[i])
            jj = vertices.index(h[j])
            graph.add_edge(ii, jj, d)
            graph.add_edge(jj, ii, d)

    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            vis = True
            # Check to see if (i, j) intersects with (k1, k2) or not, along with other conditions for visibility
            for poly in all_polys:

                for k in range(len(poly)):
                    k1 = k
                    k2 = (k + 1) % (len(poly))

                    edge = [poly[k1], poly[k2]]
                    if not is_actual_edge(all_polys, edge):
                        continue

                    # Intersect with some edges
                    if intersectProp(vertices[i], vertices[j], poly[k1], poly[k2]):
                        vis = False
                        break

                if not vis:
                    break

            edge_ij = [vertices[i], vertices[j]]

            # Exclude the exterior diagonals of the main polygon
            if vertices[i] in main_poly and vertices[j] in main_poly:
                # Not inter diagonal
                if not diagonal(i, j, main_poly):
                    continue

            # Exclude the interior diagonals of the holes
            if vertices[i] not in main_poly and is_same_poly(all_polys, edge_ij):
                # Not external diagonal
                is_diag = False
                for h in holes:
                    if vertices[i] in h and vertices[j] in h:
                        if diagonal(h.index(vertices[i]), h.index(vertices[j]), h):
                            is_diag = True
                            break
                if is_diag:
                    continue

            if vis:
                d = l2_dist(vertices[i], vertices[j])
                graph.add_edge(i, j, d)
                graph.add_edge(j, i, d)

    return graph


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


class PolyNode:
    """
    This is the class for representing each convex component as a node
    Any two components that share an edge are supposed to also share an edge in this graph
    """

    def __init__(self, poly: List):
        self.verts = copy.deepcopy(poly)
        self.edge = []
        self.visited = False

    def length(self):
        l = 0
        for i in range(len(self.verts)):
            l += l2_dist(self.verts[i], self.verts[(i + 1) % len(self.verts)])
        return l


def share_edge(poly1: PolyNode, poly2: PolyNode):
    """
    Check if two poly nodes share an edge
    :param poly1: PolyNode
    :param poly2: PolyNode
    :return: Boolean, True if they share one, False otherwise
    """
    for i in range(len(poly1.verts)):
        for j in range(len(poly2.verts)):
            i1 = i
            i2 = (i + 1) % len(poly1.verts)
            j1 = j
            j2 = (j + 1) % len(poly2.verts)
            if poly1.verts[i1] == poly2.verts[j1] and poly1.verts[i2] == poly2.verts[j2]:
                return True
            if poly1.verts[i1] == poly2.verts[j2] and poly1.verts[i2] == poly2.verts[j1]:
                return True

    return False


def create_edge(poly_nodes: List[PolyNode]):
    """
    Check if two polygons share an edge.
    If they do, create an undirected edge (two directed edges)
    :return: None
    """
    for i in range(len(poly_nodes)):
        for j in range(i, len(poly_nodes)):
            if i == j:
                continue
            if share_edge(poly_nodes[i], poly_nodes[j]):
                poly_nodes[i].edge.append(poly_nodes[j])
                poly_nodes[j].edge.append(poly_nodes[i])


def poly_decomp_cgal(verts: List[List[float]]) -> List[List[List[float]]]:
    """
    Decompose a polygon into convex components (optimally) using CGAL, coordinates are int.
    It calls an executable compile from C++.
    The path to the executable should be in poly_decomp_path global variable
    :param verts: List[List[float]], list of vertices of the polygon, must be in counterclockwise order
    :return: List[List[List[float]]], list of convex subpolygons
    """
    polygons = []
    arg = poly_decomp_path + ' --verts='
    for v in verts:
        arg += str(v[0]) + ',' + str(v[1]) + ','
    arg = arg[:-1]
    print('Running ' + arg)

    if os.name == 'posix':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True)
    if os.name == 'nt':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE)

    popen.wait()
    output = popen.stdout.read().decode("utf-8")
    output = output.split()
    output = [int(o) for o in output]

    n_poly = output[0]
    output = output[1:]
    for i in range(n_poly):
        n_vert = output[0]
        output = output[1:]
        polygon = []
        for j in range(n_vert):
            polygon += [[output[j * 2], output[j * 2 + 1]]]
        output = output[n_vert * 2:]
        polygons += [polygon]
    return polygons


def polygon_triangulation_cgal(polygon: List[List[float]], holes: List[List[List[float]]]=None) \
        -> List[List[List[float]]]:
    """
    Triangulate a polygon (possibly with holes) using CGAL. It calls an executable compile from C++.
    The path to the executable should be in poly_triangulation_path global variable
    :param polygon: List[List[float]], the polygon
    :param holes:List[List[List[float]]], list of holes
    :return: List[List[List[float]]], list of triangles
    """
    triangles = []
    arg = poly_triangulation_path + ' --verts='
    arg += str(len(polygon)) + ','
    for v in polygon:
        arg += str(v[0]) + ',' + str(v[1]) + ','

    for h in holes:
        arg += str(len(h)) + ','
        for v in h:
            arg += str(v[0]) + ',' + str(v[1]) + ','

    arg = arg[:-1]
    print('Running ' + arg)

    if os.name == 'posix':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True)
    if os.name == 'nt':
        popen = subprocess.Popen(arg, stdout=subprocess.PIPE)

    popen.wait()
    output = popen.stdout.read().decode("utf-8")
    output = output.split(',')

    for o in output:
        triangle_s = o.split()
        triangle = []
        for i in range(int(len(triangle_s) / 2)):
            triangle.append([int(triangle_s[i * 2]), int(triangle_s[i * 2 + 1])])
        if len(triangle) == 0:
            continue
        triangles.append(triangle)

    return triangles


def find_root_node(poly_nodes: List[PolyNode], agents: List[List[float]]) -> PolyNode:
    """
    Find the subpolygon in which the two agents are
    :param poly_nodes: list of all polygonal nodes
    :param agents: positions of the two agents (an edge)
    :return: PolyNode, root node
    """
    for node in poly_nodes:
        for i in range(len(node.verts)):
            i1 = i
            i2 = (i + 1) % len(node.verts)
            if (node.verts[i1] == agents[0] and node.verts[i2] == agents[1]):
                return node
            if (node.verts[i1] == agents[1] and node.verts[i2] == agents[0]):
                return node

    return PolyNode([])


def edge_in_poly(poly: PolyNode, edge: List[List[float]]):
    """
    Check if an edge (two 2D points) belongs to a PolyNode
    :param poly: PolyNode
    :param edge: List of two 2D points
    :return: True if yes, otherwise False
    """
    for i in range(len(poly.verts)):
        i1 = i
        i2 = (i + 1) % len(poly.verts)
        if poly.verts[i1] == edge[0] and poly.verts[i2] == edge[1]:
            return True
        if poly.verts[i1] == edge[1] and poly.verts[i2] == edge[0]:
            return True
    return False


def append_correct(schedule: List[List[float]], new_edge: List[float]):
    """
    Append a new edge to the schedule correctly.
    This solves the case in which two agents have to cross themselves to move on to the new edge (they shouldn't)
    :param schedule: old schedule (i.e. list of edges)
    :param new_edge: new edge, a pair of 2 floats
    :return:
    """
    old_edge = schedule[-1]
    if intersectProp(old_edge[0], new_edge[0], old_edge[1], new_edge[1]):
        new_edge[0], new_edge[1] = new_edge[1], new_edge[0]
    elif old_edge[0] == new_edge[1] or old_edge[1] == new_edge[0]:
        new_edge[0], new_edge[1] = new_edge[1], new_edge[0]
    schedule.append(new_edge)


def compute_sweep_time(i1, j1, i2, j2, verts):
    """
    Compute sweeping time from edge (i1, j1) to edge (i2, j2), in the convex polygon depicted in verts
    i1 < j1 < j2 < i2 < i1 (counter clockwise):
    i1---...---i2
    |           |
    |           |
    j1---...---j2
    :param i1: integer index
    :param j1: integer index
    :param i2: integer index
    :param j2: integer index
    :param verts: list of the polygon vertices
    :return: float, sweeping time
    """
    d1 = 0
    n = len(verts)
    while j1 != j2:
        d1 += l2_dist(verts[j1], verts[(j1 + 1) % n])
        j1 = (j1 + 1) % n

    d2 = 0
    while i2 != i1:
        d2 += l2_dist(verts[i2], verts[(i2 + 1) % n])
        i2 = (i2 + 1) % n

    return max(d1, d2)


def edge_moving_distance(e1, root, neighbor):
    """
    Compute the edge moving (sweeping) distance from e1 to the incidence edge between its root and its neighbor
    :param e1: the edge
    :param root: current node
    :param neighbor: child node
    :return: float, moviing distance (sweeping time)
    """
    i1 = root.verts.index(e1[0])
    j1 = root.verts.index(e1[1])
    if (i1 + 1) % len(root.verts) != j1:
        i1, j1 = j1, i1

    for k in range(len(neighbor.verts)):
        e2 = [neighbor.verts[k], neighbor.verts[(k + 1) % len(neighbor.verts)]]
        if edge_in_poly(root, e2):
            i2 = root.verts.index(e2[0])
            j2 = root.verts.index(e2[1])
            if (i2 - 1) % len(root.verts) != j2:
                i2, j2 = j2, i2

            return compute_sweep_time(i1, j1, i2, j2, root.verts)


def find_cloest_neighbor(root, e, neighbors):
    """
    Find the neighbor of root that has the smallest sweeping time from edge e
    :param root: current polygon
    :param e: the edge to query
    :param neighbors: all neighbors of root
    :return: int, index of the best neighbor
    """
    if len(neighbors) == 1:
        return 0

    min_ind = 0
    min_dist = edge_moving_distance(e, root, neighbors[0])
    for i in range(1, len(neighbors)):
        new_dist = edge_moving_distance(e, root, neighbors[i])
        if new_dist < min_dist:
            min_dist = new_dist
            min_ind = i

    return min_ind


def dfs(polygons: List[PolyNode], root: PolyNode,
        polygon: List[List[float]], vis_graph: Graph,
        agents: List[List[float]], unvisited_children: List[PolyNode]) -> List[List[List[List[float]]]]:
    """
    The main function for sweeping with two agents.
    :param polygons: List[PolyNode], list of PolyNode in which each is a subpolygon after partitioning
    :param root: PolyNode, first subpolygon that the two agents start in
    :param polygon: List[List[float]], original polygon
    :param vis_graph: Graph, the visibility graph, used to find shortest path from a node to another
    :param agents: List[List[float]], starting locations of the two agents
    :param unvisited_children: List[PolyNode], list of unvisited branches so far in the DFS traversal
    :return: List[List[List[List[float]]]] (4 List), the strategy. The list of states of the two agents (line segments)
    Each of the inner List[List[List[float]]] (3 List) is a completely valid strategy, consisting of multiple line segments.
    """

    if root.visited:
        return []

    if root in unvisited_children:
        del unvisited_children[unvisited_children.index(root)]

    root.visited = True
    # print(str(root.verts) + ' ' + str(polygons.index(root)))
    schedule = []

    # Generate a list of children nodes, in counterclockwise order
    # edges is the list of edges, same length as children list
    # Each edge here is the edge that connects current node with the corresponding child
    children = []
    edges = []
    for i in range(len(root.verts)):
        i1 = i
        i2 = (i + 1) % len(root.verts)
        edge = [root.verts[i1], root.verts[i2]]
        for poly in polygons:
            if poly == root or poly.visited:
                continue
            if edge_in_poly(poly, edge):
                children.append(poly)
                edges.append(edge)

    # Generate a sweep schedule for the current root node (a convex subpoly)
    # Note that at the end of the schedule, the two agents will align with the edge shared with the best child found in find_closest_neighbor

    # Find shortest path from the previous locations of the agents to the gate of this node
    if (agents[0] not in root.verts) or (agents[1] not in root.verts):
        best_d = float('inf')
        best_convergent_point = -1
        best_next_point = -1

        leaf_poly = None
        for poly in polygons:
            if edge_in_poly(poly, [agents[0], agents[1]]):
                leaf_poly = poly
                break

        for e in range(len(root.verts)):
            next_edge = [polygon.index(root.verts[e]), polygon.index(root.verts[(e + 1) % len(root.verts)])]

            min_d = float('inf')
            convergent_point = -1
            next_point = -1

            for i in range(len(leaf_poly.verts)):
                v = leaf_poly.verts[i]

                dijkstra = DijkstraSPF(vis_graph, polygon.index(v))
                path_lengths = [dijkstra.get_distance(next_edge[0]), dijkstra.get_distance(next_edge[1])]
                total_length = min(path_lengths)
                total_length += max(l2_dist(v, agents[0]), l2_dist(v, agents[1]))
                if total_length < min_d:
                    min_d = total_length
                    convergent_point = polygon.index(v)

                    if path_lengths[0] < path_lengths[1]:
                        next_point = next_edge[0]
                    else:
                        next_point = next_edge[1]

            if min_d < best_d:
                best_convergent_point = convergent_point
                best_d = min_d

                best_next_point = next_point

        schedule_r = [agents]
        schedule_r += [[polygon[best_convergent_point], polygon[best_convergent_point]]]
        schedule += [schedule_r]

        dijkstra = DijkstraSPF(vis_graph, best_convergent_point)
        path = dijkstra.get_path(best_next_point)
        schedule_r = []
        for p in path:
            schedule_r += [[polygon[p], polygon[p]]]

        schedule += [schedule_r]

        agents = schedule_r[-1]

    schedule_r = []
    # Special case if leaf node
    if len(children) == 0:
        l = l2_dist(agents[0], agents[1])
        poly_len = root.length()

        # i1 < j1 < j2 < i2 < i1 (counter clockwise)
        # i1---...---i2
        # |           |
        # |           |
        # j1---...---j2
        i1 = root.verts.index(agents[0])
        j1 = root.verts.index(agents[1])
        if (i1 + len(root.verts) - 1) % len(root.verts) == j1:
            i1, j1 = j1, i1
        schedule_r.append([root.verts[i1], root.verts[j1]])

        best_edge_ind = -1
        best_sweep_time = float('inf')
        for k in range(len(root.verts)):
            j2 = k
            i2 = (k + 1) % len(root.verts)
            if j2 == i1:
                continue
            sweep_time = compute_sweep_time(copy.copy(i1), copy.copy(j1), i2, j2, root.verts)
            if sweep_time < best_sweep_time:
                best_sweep_time = sweep_time
                best_edge_ind = k

        j2 = best_edge_ind
        i2 = (best_edge_ind + 1) % len(root.verts)
        di = 0
        dj = 0
        while i1 != i2 or j1 != j2:
            if i1 != i2:
                new_di = di + l2_dist(root.verts[(i1 - 1) % len(root.verts)], root.verts[i1])
            else:
                new_di = di
                j1 = (j1 + 1) % len(root.verts)
                schedule_r.append([root.verts[i1], root.verts[j1]])
                continue

            if j1 != j2:
                new_dj = dj + l2_dist(root.verts[(j1 + 1) % len(root.verts)], root.verts[j1])
            else:
                new_dj = dj
                i1 = (i1 - 1) % len(root.verts)
                schedule_r.append([root.verts[i1], root.verts[j1]])
                continue

            if new_di < new_dj:
                i1 = (i1 - 1) % len(root.verts)
            else:
                j1 = (j1 + 1) % len(root.verts)
            schedule_r.append([root.verts[i1], root.verts[j1]])

        # append_correct(schedule_r, [root.verts[j2], root.verts[i2]])

        schedule += [schedule_r]

        if len(unvisited_children) == 0:
            return schedule

        best_d_c = float('inf')
        best_c = -1
        for c in range(len(unvisited_children)):
            child = unvisited_children[c]
            best_d = float('inf')
            for e in range(len(child.verts)):
                next_edge = [polygon.index(child.verts[e]), polygon.index(child.verts[(e + 1) % len(child.verts)])]

                min_d = float('inf')

                for i in range(len(root.verts)):
                    v = root.verts[i]

                    dijkstra = DijkstraSPF(vis_graph, polygon.index(v))
                    path_lengths = [dijkstra.get_distance(next_edge[0]), dijkstra.get_distance(next_edge[1])]
                    total_length = min(path_lengths)
                    total_length += max(l2_dist(v, agents[0]), l2_dist(v, agents[1]))
                    if total_length < min_d:
                        min_d = total_length
                        convergent_point = polygon.index(v)

                        # if path_lengths[0] < path_lengths[1]:
                        #     next_point = next_edge[0]
                        # else:
                        #     next_point = next_edge[1]

                if min_d < best_d:
                    # best_convergent_point = convergent_point
                    best_d = min_d

                    # best_next_point = next_point
            if best_d < best_d_c:
                best_d_c = best_d
                best_c = c

        schedule_r = dfs(polygons, unvisited_children[best_c], polygon, vis_graph, schedule_r[-1], unvisited_children)
        # del unvisited_children[best_c]
        if len(schedule_r) > 1:
            schedule += schedule_r

        return schedule

    # General case, moving the segment from edge/diag (i1, j1) to edge/diag (i2, j2)
    i1 = root.verts.index(agents[0])
    j1 = root.verts.index(agents[1])

    closest_neighbor = find_cloest_neighbor(root, agents, children)

    if (i1 + len(root.verts) - 1) % len(root.verts) == j1:
        i1, j1 = j1, i1

    i2 = root.verts.index(edges[closest_neighbor][0])
    j2 = root.verts.index(edges[closest_neighbor][1])
    if (j2 + len(root.verts) - 1) % len(root.verts) == i2:
        i2, j2 = j2, i2

    schedule_r.append([root.verts[i1], root.verts[j1]])

    while i1 != i2 or j1 != j2:
        if i1 != i2:
            i1 = (i1 + len(root.verts) - 1) % len(root.verts)
        if j1 != j2:
            j1 = (j1 + 1) % len(root.verts)

        schedule_r.append([root.verts[i1], root.verts[j1]])

    schedule += [schedule_r]

    last_edge = edges[closest_neighbor]

    while len(children) > 0:
        # Find the next child with start_edge of closest moving distance from last_edge
        c = children[closest_neighbor]
        del children[closest_neighbor]
        closest_neighbor = 0
        if not c.visited:
            unvisited_children += children
            schedule_c = dfs(polygons, c, polygon, vis_graph, last_edge, unvisited_children)

            # del edges[closest_neighbor]
            if len(schedule_c) > 0:
                last_edge = schedule_c[-1][-1]
                schedule += schedule_c

    return schedule


def draw_polys(ps, ax, linewidth=0.5, edgecolor='0.7'):
    """
    Draw polygons
    :param ps: List[PolyNode]
    :param ax: ax object created by matplotlib
    :param linewidth: float, width of the lines used to draw the polygons
    :param edgecolor: float (grayscale), or [float, float, float] (RGB), color for the polygon edges
    :return: None
    """
    for p in ps:
        p_np = np.asarray(p.verts)
        poly = patches.Polygon(p_np, fill=False, edgecolor=edgecolor, linewidth=linewidth)
        ax.add_patch(poly)
        # for v in p.verts:
        #   plt.text(v[0], v[1]+1, str(v[0]) + ',' + str(v[1]))

    for p in ps:
        p_np = np.asarray(p.verts)
        plt.scatter(p_np[::-1, 0], p_np[::-1, 1])

    plt.draw()


def draw_vis_graph(vertices: List[List[float]], g: Graph) -> None:
    """
    Draw the visibility graph of the vertices in a polygon
    :param vertices: List[List[float]]], list of vertices
    :param g: Graph, the graph computed from create_vis_graph or create_vis_graph_with_holes
    :return: None
    """
    for node in list(g.get_nodes()):
        for neighbor in list(g.get_adjacent_nodes(node)):
            v1 = vertices[node]
            v2 = vertices[neighbor]
            plt.plot([v1[0], v2[0]], [v1[1], v2[1]], color='r', linewidth=0.5)
            if node != 30:
                plt.text(v1[0] + 0.5, v1[1], str(node))
            if neighbor != 30:
                plt.text(v2[0] + 0.5, v2[1], str(neighbor))


def get_curr_point(verts, t, speed):
    """
    verts store two points a and b, an agent will move from a to b at speed 'speed'.
    This function computes the location of the agent at time t.
    :param verts: List[List[float]]
    :param t: float, time
    :param speed: float, speed
    :return: List[float], location of the agent
    """
    total_time = 0
    for i in range(len(verts) - 1):
        p1 = np.array(verts[i]).astype(float)
        p2 = np.array(verts[i + 1]).astype(float)
        total_time_ = total_time + l2_dist(verts[i], verts[i + 1]) / speed
        if total_time_ > t:
            ratio = (t - total_time) / (total_time_ - total_time)
            new_p = (ratio * (p2 - p1)) + p1
            return new_p

        total_time = total_time_

    return verts[-1]


def animate_schedule(polygons, subpolygons, schedule, speed):
    """
    Animate the computed sweeping strategy
    :param polygon: List[List[float]], the original polygon and its holes
    :param polygons: List[PolyNode], list of convex subpolygons
    :param schedule: List[List[List[List[float]]]], the strategy.
    :param speed: float, animation speed
    :return: None
    """
    fig, ax = plt.subplots()
    draw_polys([PolyNode(p) for p in polygons], ax, edgecolor=[0, 0, 0], linewidth=1)
    draw_polys(subpolygons, ax)
    plt.axis('equal')
    plt.waitforbuttonpress()

    dt = 0.001
    ela_time = 0
    a = schedule[0][0][0]
    b = schedule[0][0][1]
    line = ax.plot([a[0], b[0]], [a[1], b[1]], color='b')[0]
    scatter = ax.scatter([a[0], b[0]], [a[1], b[1]], edgecolors='b')

    for s in schedule:
        t0 = time.time()
        path_a = [x[0] for x in s]
        path_b = [x[1] for x in s]
        while True:
            t1 = time.time()
            ela_time = t1 - t0
            a = get_curr_point(path_a, ela_time, speed)
            b = get_curr_point(path_b, ela_time, speed)

            line.set_data([a[0], b[0]], [a[1], b[1]])
            scatter.set_offsets([[a[0], a[1]], [b[0], b[1]]])
            plt.show()
            plt.axis('equal')
            plt.pause(0.01)

            try:
                if list(a) == s[-1][0] and list(b) == s[-1][1]:
                    break
            except Exception as e:
                abc = 1


def animate_schedule_gif(polygons, subpolygons, schedule, speed, boundary_polygon):
    """
    Animate the computed sweeping strategy
    :param polygon: List[List[float]], the original polygon and its holes
    :param polygons: List[PolyNode], list of convex subpolygons
    :param schedule: List[List[List[List[float]]]], the strategy.
    :param speed: float, animation speed
    :return: None
    """
    fig, ax = plt.subplots()
    plt.plot(boundary_polygon[:, 0], boundary_polygon[:, 1], 'k')
    rectangle = np.min(boundary_polygon, 0)
    r_width, r_height = np.max(boundary_polygon, 0) - np.min(boundary_polygon, 0)
    patch = patches.Rectangle(rectangle, r_width, r_height, facecolor = [1,1,1])
    ax.add_patch(patch)
    #draw_polys([PolyNode(p) for p in polygons], ax, edgecolor=[0, 0, 0], linewidth=1)
    #draw_polys(subpolygons, ax)
    plt.axis('equal')
    plt.waitforbuttonpress()

    dt = 0.001
    ela_time = 0
    a = schedule[0][0][0]
    b = schedule[0][0][1]
    line = ax.plot([a[0], b[0]], [a[1], b[1]], color='b')[0]
    scatter = ax.scatter([a[0], b[0]], [a[1], b[1]], edgecolors='b', zorder=100)

    segment_path = np.asarray([[[a[0], b[0]], [a[1], b[1]]]])
    #print(segment_path.shape)
    segment_path = np.concatenate((segment_path, segment_path), 0)
    for s in schedule:
        t0 = time.time()
        path_a = [x[0] for x in s]
        path_b = [x[1] for x in s]
        while True:
            t1 = time.time()
            ela_time = t1 - t0
            a = get_curr_point(path_a, ela_time, speed)
            b = get_curr_point(path_b, ela_time, speed)
            segment_path[1, :, :] = np.asarray([[a[0], b[0]], [a[1], b[1]]])

            x = [segment_path[0, 0, 0], segment_path[0, 0, 1], segment_path[1, 0, 1], segment_path[1, 0, 0]]
            y = [segment_path[0, 1, 0], segment_path[0, 1, 1], segment_path[1, 1, 1], segment_path[1, 1, 0]]
            plt.fill(x, y, facecolor = [1, 0.8, 0.5], alpha=1)
            line.set_data([a[0], b[0]], [a[1], b[1]])
            scatter.set_offsets([[a[0], a[1]], [b[0], b[1]]])
            plt.show()
            plt.axis('equal')
            plt.pause(0.001)

            segment_path[0, :, :] = segment_path[1, :, :]

            try:
                if list(a) == s[-1][0] and list(b) == s[-1][1]:
                    break
            except Exception as e:
                abc = 1


def test_without_hole():
    plt.ion()
    speed = 15

    # polygon = [[0, 0], [10, 0], [10, 10], [20, 20], [17, 28], [15, 30], [7, 15], [3, 15], [-5, 30], [-10, 20], [0, 10]]
    # agents = [[0, 0], [10, 0]]

    # polygon = [[0, 0], [37, 0], [37, 8], [36, 11], [36, 22], [33, 22], [28, 24], [27, 20], [31, 21],
    #           [30, 19], [20, 19], [22, 24], [14, 16], [20, 14], [30, 14],
    #           [32, 11], [0, 10]]
    # agents = [[0, 0], [0, 10]]

    # polygon = [[0, 0], [40, 0], [40, 10], [38, 10], [38, 20]]
    # agents = [[0, 10], [0, 0]]

    # speed = 40
    # polygon = [[40, 0], [0, 60], [10, 50], [22, 55], [20, 45], [30, 42], [27, 50],
    #            [40, 40], [53, 50], [50, 42], [60, 45], [58, 55], [70, 50], [80, 60]]
    # polygon = polygon[::-1]
    # agents = [[10, 50], [0, 60]]

    # Random 1
    # speed = 10
    # polygon = [[0, 0], [1, 3], [2, 1], [3, 5], [6, 8], [13, 4], [13, 0], [14, -1],[14, 3],
    #            [20, 2], [17, 4], [21, 6], [14, 5], [7, 9], [9, 10], [8, 16],
    #            [13, 15], [14, 17], [18, 17], [16, 18], [17, 19], [8, 18], [4, 18], [6, 15], [4, 15],
    #            [5, 10], [-5, 7], [1, 8], [1, 5]]
    # #polygon = polygon[::-1]
    # agents = [[0, 0], [1, 3]]

    # Flower
    # speed = 10
    # polygon = [[0, 0], [4, 6], [8, 0], [6, 7], [13, 5], [7, 9],
    #            [13, 13], [6, 11], [8, 18], [4, 12], [0, 18],
    #            [2, 11], [-5, 13], [1, 9], [-5, 5], [2, 7]
    #            ]
    # agents = [[0,0], [4, 6]]

    # pinwheel
    # speed = 10
    # polygon = [[3, -1], [0, 7], [10, 0], [2, 6], [13, 10], [4, 7],
    #            [9, 15], [5, 9], [1, 19], [4, 11], [-6, 17],
    #            [2, 12], [-8, 8], [0, 11], [-5, 3], [-1, 9]
    #            ]
    # agents = [[3, -1], [0, 7]]

    # random 2
    # speed = 150
    # polygon = [
    #     [480, 720], [352, 576], [304, 592], [272, 528], [288, 752], [32, 720],
    #     [48, 688], [144, 624], [192, 704], [256, 704], [224, 528], [64, 496],
    #     [80, 432], [336, 480], [336, 368], [448, 352], [368, 496], [416, 608],
    #     [496, 576], [448, 272], [544, 272], [512, 496], [576, 464], [496, 624],
    #     [560, 656], [560, 800], [272, 768], [368, 736], [528, 784], [464, 656]
    # ]
    # agents = [[480, 720], [352, 576]]

    # polygon = [
    #     [0,0], [10,0], [11,-1],[12,0], [22,0],
    #     [22,10], [12,10], [11,11], [10, 10], [0, 10]
    # ]
    # agents = [[0, 0], [0, 10]]

    # Random 1
    speed = 200
    polygon = [[240, 672], [128, 592], [208, 544], [80, 496],
               [96, 400], [176, 336], [192, 512], [304, 448],
               [304, 320], [432, 272], [592, 384], [464, 496],
               [368, 384], [352, 592], [528, 608], [432, 784],
               [112, 800], [96, 752], [384, 704], [256, 512]]
    agents = [[240, 672], [128, 592]]

    # Parition the polygon into conve pieces
    polygons = poly_decomp_cgal(polygon)
    # Create visibility graph between vertices
    vis_graph = create_vis_graph(polygon, [])

    polygons = [PolyNode(poly) for poly in polygons]

    # Create a graph in which each subpolygon is a node, two nodes share an edge if two corresponding subpolygons share
    # an edge in the partition
    create_edge(polygons)
    # Find the first subpolygon that the agents reside in
    root = find_root_node(polygons, agents)
    # Main algorithm
    schedule = dfs(polygons, root, polygon, vis_graph, agents, [])

    animate_schedule([polygon], polygons, schedule, speed)
    plt.waitforbuttonpress()
    plt.close()


def test_with_hole():
    plt.ion()

    # speed = 7
    # polygon = [[0,0], [10, 0], [10, 10], [0, 10]]
    # agents = [[0, 0], [0, 10]]
    # holes = [[[1, 1], [3, 1], [3, 3], [1, 3]],
    #          [[4, 4], [6, 4], [6, 6], [5, 5], [4, 6]]]

    # Random 1
    speed = 100
    polygon = [[240, 672], [128, 592], [208, 544], [80, 496],
               [96, 400], [176, 336], [192, 512], [304, 448],
               [304, 320], [432, 272], [592, 384], [464, 496],
               [368, 384], [352, 592], [528, 608], [432, 784],
               [112, 800], [96, 752], [384, 704], [256, 512]]
    agents = [[240, 672], [128, 592]]
    holes = [
        [[112, 496], [96, 464], [160, 400], [176, 464],
         [160, 480], [144, 432], [128, 448], [144, 480]],  # hole 1
        [[288, 496], [336, 432], [320, 352], [416, 304],
         [464, 416], [544, 368], [528, 416], [432, 432],
         [400, 368], [352, 368], [336, 480]],  # hole 2
        [[368, 768], [416, 656], [368, 640], [496, 624],
         [480, 672], [432, 640], [450, 704], [448, 672],
         [464, 704], [416, 720], [432, 736], [390, 736]],  # hole 3
        [[160, 784], [352, 768], [352, 720], [320, 752],
         [208, 752]]  # hole 4
    ]
    holes[3] = holes[3][::-1]

    # speed = 200
    # polygon = [[64, 704], [64, 320], [512, 320], [512, 704]]
    # agents = [[64, 704], [64, 320]]
    # holes = [
    #     [[126, 640], [126, 560], [193, 558], [194, 466],
    #      [126, 464], [126, 384], [208, 384], [210, 446],
    #      [366, 450], [368, 384], [450, 384], [450, 464],
    #      [386, 466], [384, 558], [452, 562], [450, 640],
    #      [368, 640], [366, 578], [210, 578], [208, 640]],       # hole 0
    #     [[110, 544], [170, 544], [170, 480], [110, 480]],       # hole 1
    #     [[222, 352], [222, 432], [354, 432], [354, 352],
    #      [334, 352], [334, 416], [242, 416], [242, 352]],       # hole 2
    #     [[258, 400], [322, 400], [322, 368], [258, 368]],       # hole 3
    #     [[402, 544], [402, 480], [467, 480], [467, 544]],       # hole 4
    #     [[240, 592], [240, 660], [352, 660], [352, 592],
    #      [336, 592], [336, 640], [256, 640], [256, 592]],       # hole 5
    #     [[272, 624], [272, 608], [320, 608], [320, 624]]        # hole 6
    # ]
    # holes[1] = holes[1][::-1]
    # holes[2] = holes[2][::-1]
    # holes[3] = holes[3][::-1]
    # #holes[4] = holes[4][::-1]
    # holes[5] = holes[5][::-1]

    # Triangulation
    triangles = polygon_triangulation_cgal(polygon, holes)
    # Compute visibility graph between vertices graph
    vis_graph = create_vis_graph_with_holes(polygon, holes)

    # This is to visualize the vis graph to see if it's correct
    # vertices = copy.deepcopy(polygon)
    # for h in holes:
    #     vertices += h
    # draw_vis_graph(vertices, vis_graph)

    triangles = [PolyNode(t) for t in triangles]

    # Create a graph in which each triangle is a node, two nodes share an edge if their triangles share an edge
    create_edge(triangles)
    # Find the first node that the two agents reside in
    root = find_root_node(triangles, agents)

    all_vert = copy.deepcopy(polygon)
    for h in holes:
        all_vert += h
    # Main algorithm
    schedule = dfs(triangles, root, all_vert, vis_graph, agents, [])
    animate_schedule([polygon] + holes, triangles, schedule, speed)
    plt.waitforbuttonpress()
    plt.close()


def create_gif():
    plt.ion()

    # Problem definition gif
    # speed = 100
    # polygon = [[128, 608], [128, 612], [160, 640], [224, 608], [272, 656],
    #            [336, 624], [304, 556], [380, 568], [336, 500], [256, 516],
    #            [260, 588], [224, 576], [160, 608]]
    # polygon = polygon[::-1]
    # boundary_polygon = [[108, 664], [128, 556], [204, 560], [240, 500], [368, 484],
    #                     [436, 560], [404, 604], [340, 596], [352, 656], [224, 676],
    #                     [188, 652], [108, 664]]
    # boundary_polygon = np.array(boundary_polygon)
    # agents = [[128, 608], [128, 612]]

    # CG gif
    speed = 70
    polygon = [[160, 624],[176, 624],[176, 640],[160, 656],[128, 656],
               [112, 640],[112, 576],[128, 560],[160, 560],[176, 576],
               [176, 584],[192, 584],[192, 576],[208, 560],[240, 560],
               [256, 576],[256, 592],[264, 592],[264, 600],[232, 600],
               [232, 592],[240, 592],[240, 576],[208, 576],[208, 640],
               [240, 640],[240, 624],[256, 624],[256, 640],[240, 656],
               [208, 656],[192, 640],[192, 588],[176, 588],[176, 592],
               [160, 592],[160, 576], [128, 576],[128, 640], [160, 640]]
    #polygon = polygon[::-1]
    boundary_polygon = [[288, 576], [288, 640], [256, 672], [128, 672], [96, 640], [96, 576], [128, 544], [256, 544], [288, 576]]
    boundary_polygon = np.array(polygon)
    boundary_polygon = np.concatenate((boundary_polygon, np.asarray([[160, 624]])), 0)

    agents = [[160, 624], [160, 640]]

    # Parition the polygon into convex pieces
    polygons = poly_decomp_cgal(polygon)

    # Create visibility graph between vertices
    vis_graph = create_vis_graph(polygon, [])

    polygons = [PolyNode(poly) for poly in polygons]

    # Create a graph in which each subpolygon is a node, two nodes share an edge if two corresponding subpolygons share
    # an edge in the partition
    create_edge(polygons)
    # Find the first subpolygon that the agents reside in
    root = find_root_node(polygons, agents)
    # Main algorithm
    schedule = dfs(polygons, root, polygon, vis_graph, agents, [])
    animate_schedule_gif([polygon], polygons, schedule, speed, boundary_polygon)
    plt.waitforbuttonpress()
    plt.close()


if __name__ == '__main__':
    create_gif()
    #test_without_hole()
    #test_with_hole()
