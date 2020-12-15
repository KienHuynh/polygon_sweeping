import subprocess
import poly_decomp as pd

from typing import List
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import numpy as np
from dijkstra import Graph, DijkstraSPF
import copy
import sys
import time

from base import *

poly_decomp_path = 'D:/study/phd_research/main/polygon_sweeping/poly_decomp_cgal/Debug/poly_decomp_cgal.exe'


# Will need to change this later
# class Graph():
#     def __init__(self, vertices):
#         self.V = range(len(vertices))
#         self.graph = [[0.0] * len(vertices) for x in vertices]
#
#
#
#     def stp(self):
#         abc = 1

def create_vis_graph(vertices):
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


def l2_dist(p1, p2):
    p1 = np.asarray(p1).astype(float)
    p2 = np.asarray(p2).astype(float)
    return float(np.sqrt(np.sum((p1 - p2)**2)))


class PolyNode:
    def __init__(self, poly: List):
        self.verts = copy.deepcopy(poly)
        self.edge = []
        self.visited = False

    def length(self):
        l = 0
        for i in range(len(self.verts)):
            l += l2_dist(self.verts[i], self.verts[(i+1) % len(self.verts)])
        return l



def share_edge(poly1: PolyNode, poly2: PolyNode):
    for i in range(len(poly1.verts)):
        for j in range(len(poly2.verts)):
            i1 = i
            i2 = (i+1) % len(poly1.verts)
            j1 = j
            j2 = (j+1) % len(poly2.verts)
            if poly1.verts[i1] == poly2.verts[j1] and poly1.verts[i2] == poly2.verts[j2]:
                return True
            if poly1.verts[i1] == poly2.verts[j2] and poly1.verts[i2] == poly2.verts[j1]:
                return True

    return False


def create_edge(poly_nodes: List[PolyNode]):
    """
    Check if two polygons share an edge.
    If they do, create an undirected edge (two directed edges)
    :return:
    """
    for i in range(len(poly_nodes)):
        for j in range(i, len(poly_nodes)):
            if i == j:
                continue
            if share_edge(poly_nodes[i], poly_nodes[j]):
                poly_nodes[i].edge.append(poly_nodes[j])
                poly_nodes[j].edge.append(poly_nodes[i])


def poly_decomp_cgal(verts):
    """ Decompose a polygon into convex components (optimally), coordinates are int
    """
    polygons = []
    arg = poly_decomp_path + ' --verts='
    for v in verts:
        arg += str(v[0]) + ',' + str(v[1]) + ','
    arg = arg[:-1]
    print('Running ' + arg)

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
            polygon += [[output[j*2], output[j*2+1]]]
        output = output[n_vert*2:]
        polygons += [polygon]
    return polygons


def find_root_node(poly_nodes: List[PolyNode], agents: List[List[float]]) -> PolyNode:
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
    for i in range(len(poly.verts)):
        i1 = i
        i2 = (i+1) % len(poly.verts)
        if poly.verts[i1] == edge[0] and poly.verts[i2] == edge[1]:
            return True
        if poly.verts[i1] == edge[1] and poly.verts[i2] == edge[0]:
            return True
    return False


def append_correct(schedule: List[List[float]], new_edge: List[float]):
    """
    Append a new edge to the schedule correctly.
    This solves the case in which two agents have to cross themselves to move on to the new edge (they don't have to)
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


# i1 < j1 < j2 < i2 < i1 (counter clockwise):
# i1--...---i2
# |          |
# |          |
# j1--...---j2
def compute_sweep_time(i1, j1, i2, j2, verts):
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


# i1 < j1 < j2 < i2 < i1 (counter clockwise):
# i1--...---i2
# |          |
# |          |
# j1--...---j2
def edge_moving_distance(e1, root, neighbor):
    i1 = root.verts.index(e1[0])
    j1 = root.verts.index(e1[1])
    if (i1 + 1) % len(root.verts) != j1:
        i1, j1 = j1, i1

    for k in range(len(neighbor.verts)):
        e2 = [neighbor.verts[k], neighbor.verts[(k+1) % len(neighbor.verts)]]
        if edge_in_poly(root, e2):
            i2 = root.verts.index(e2[0])
            j2 = root.verts.index(e2[1])
            if (i2 - 1) % len(root.verts) != j2:
                i2, j2 = j2, i2

            # di = 0
            # if (i2 != i1):
            #     di += l2_dist(root.verts[i2], root.verts[(i2 - 1) % len(root.verts)])
            #     while (i2 - 1) % len(root.verts) != i1:
            #         i2 = (i2 - 1) % len(root.verts)
            #         di += l2_dist(root.verts[i2], root.verts[(i2 - 1) % len(root.verts)])
            #
            # dj = 0
            # if (j2 != j1):
            #     dj = l2_dist(root.verts[j2], root.verts[(j2 + 1) % len(root.verts)])
            #     while (j2 + 1) % len(root.verts) != j1:
            #         j2 = (j2 + 1) % len(root.verts)
            #         dj += l2_dist(root.verts[j2], root.verts[(j2 + 1) % len(root.verts)])

            return compute_sweep_time(i1, j1, i2, j2, root.verts)


def find_cloest_neighbor(root, e, neighbors):
    if (len(neighbors) == 1):
        return 0

    min_ind = 0
    min_dist = edge_moving_distance(e, root, neighbors[0])
    for i in range(1, len(neighbors)):
        new_dist = edge_moving_distance(e, root, neighbors[i])
        if new_dist < min_dist:
            min_dist = new_dist
            min_ind = i

    return min_ind


def dfs(polygons, root:PolyNode, polygon, vis_graph, pred_sib_edges, agents):
    root.visited = True
    schedule_r = []

    # Generate a list of children nodes, in counterclockwise order
    # edges is the list of edges, same length as children list
    # Each edge here is the edge that connects current node with the corresponding child
    children = []
    edges = []
    for i in range(len(root.verts)):
        i1 = i
        i2 = (i+1) % len(root.verts)
        edge = [root.verts[i1], root.verts[i2]]
        for poly in polygons:
            if poly == root or poly.visited:
                continue
            if edge_in_poly(poly, edge):
                children.append(poly)
                edges.append(edge)

    # Generate a sweep schedule for the current root node (a convex subpoly)
    # Note that at the end of the schedule, the two agents will align with the edge shared with the first child found in the above loop
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

        # to_cover_len = (poly_len - l) / 2.0
        # l0 = 0
        # while l0 <= to_cover_len:
        #     i2 = (i1 + len(root.verts) - 1) % len(root.verts)
        #     j2 = (j1 + 1) % len(root.verts)
        #     if (l0 + l2_dist(root.verts[j1], root.verts[j2])) > to_cover_len:
        #
        #         #schedule_r.append([root.verts[j1], root.verts[j2]])
        #         append_correct(schedule_r, [root.verts[j1], root.verts[j2]])
        #         break
        #     else:
        #         l0 = l0 + l2_dist(root.verts[j1], root.verts[j2])
        #
        #     if i2 == j2:
        #         dj = l2_dist(root.verts[j1], root.verts[j2])
        #         di = l2_dist(root.verts[i1], root.verts[i2])
        #         if di < dj:
        #             append_correct(schedule_r, [root.verts[i2], root.verts[j1]])
        #         else:
        #             append_correct(schedule_r, [root.verts[i1], root.verts[j2]])
        #         break
        #
        #     i1 = i2
        #     j1 = j2

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

        #append_correct(schedule_r, [root.verts[j2], root.verts[i2]])

        schedule = [schedule_r]
        if len(pred_sib_edges) == 0:
            return schedule

        # Find its way to the next branch
        # This is not a good way to do it
        next_edge = pred_sib_edges[0]
        next_edge = [polygon.index(e) for e in next_edge]

        min_d = float('inf')
        convergent_point = -1
        next_point = -1
        for i in range(len(root.verts)):
            v = root.verts[i]

            dijkstra = DijkstraSPF(vis_graph, polygon.index(v))
            path_lengths = [dijkstra.get_distance(next_edge[0]), dijkstra.get_distance(next_edge[1])]
            total_length = min(path_lengths)
            total_length += max(l2_dist(v, schedule_r[-1][0]), l2_dist(v, schedule_r[-1][1]))
            if total_length < min_d:
                min_d = total_length
                convergent_point = polygon.index(v)

                if path_lengths[0] < path_lengths[1]:
                    next_point = next_edge[0]
                else:
                    next_point = next_edge[1]

        schedule_r = [schedule_r[-1]] + [[polygon[convergent_point], polygon[convergent_point]]]
        schedule.append(schedule_r)
        schedule += [[[polygon[convergent_point], polygon[convergent_point]]]]

        dijkstra = DijkstraSPF(vis_graph, convergent_point)
        path = dijkstra.get_path(next_point)
        for p in path:
            schedule[-1] += [[polygon[p], polygon[p]]]
        return schedule

        # dijkstra1 = DijkstraSPF(vis_graph, polygon.index(schedule_r[-1][0]))
        # dijkstra2 = DijkstraSPF(vis_graph, polygon.index(schedule_r[-1][1]))
        # #dijkstra = [dijkstra1, dijkstra1, dijkstra2, dijkstra2]
        # #lengths = [dijkstra[0].get_distance(next_edge[0]),
        # #           dijkstra[1].get_distance(next_edge[1]),
        # #           dijkstra[2].get_distance(next_edge[0]),
        # #           dijkstra[3].get_distance(next_edge[1])]
        # #paths = [dijkstra[0].get_path(next_edge[0]),
        # #           dijkstra[1].get_path(next_edge[1]),
        # #           dijkstra[2].get_path(next_edge[0]),
        # #           dijkstra[3].get_path(next_edge[1])]
        # #dijkstra = [[dijkstra1, dijkstra1], [dijkstra2, dijkstra2]]
        # lengths = [[dijkstra1.get_distance(next_edge[0]),
        #            dijkstra1.get_distance(next_edge[1])],
        #            [dijkstra2.get_distance(next_edge[0]),
        #             dijkstra2.get_distance(next_edge[1])]
        #            ]
        # best_length = max(min(lengths[0]), min(lengths[1]))
        # best_path = []
        # if best_length == lengths[0][0]:
        #     best_path = dijkstra1.get_path(next_edge[0])
        # elif best_length == lengths[0][1]:
        #     best_path = dijkstra1.get_path(next_edge[1])
        # elif best_length == lengths[1][0]:
        #     best_path = dijkstra2.get_path(next_edge[0])
        # elif best_length == lengths[1][1]:
        #     best_path = dijkstra2.get_path(next_edge[1])
        #
        # min_ind = np.argmin(lengths)
        # if (len(best_path) != 0):
        #     path = best_path
        #     path = path[1:]
        #     path = [polygon[p] for p in path]
        #     schedule_s = [schedule_r[-1]] + [[path[0], path[0]]]
        #     schedule_geo = [[p, p] for p in path]
        #     schedule += [schedule_s, schedule_geo]
        # else:
        #     path = [polygon[next_edge[min_ind]], polygon[next_edge[min_ind]]]
        #     schedule_s = [schedule_r[-1]] + [[path[0], path[0]]]
        #     schedule += [schedule_s]
        #
        # return schedule


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

    schedule = [schedule_r]

    last_edge = edges[closest_neighbor]
    del edges[closest_neighbor]
    if len(edges) > 0:
        pred_sib_edges_ = edges
    else:
        pred_sib_edges_ = pred_sib_edges

    while len(children) > 0:
        # Find the next child with start_edge of closest moving distance from last_edge
        #closest_neighbor = find_cloest_neighbor(root, last_edge, children)
        #closest_neighbor = 0
        c = children[closest_neighbor]
        #c = children[0]

        schedule_c = dfs(polygons, c, polygon, vis_graph, pred_sib_edges_, last_edge)
        del children[closest_neighbor]
        closest_neighbor = 0
        #del edges[closest_neighbor]

        if len(edges) > 1:
            pred_sib_edges_ = edges
            del edges[closest_neighbor-1]
        else:
            pred_sib_edges_ = pred_sib_edges
            edges = edges[1:]
        last_edge = schedule_c[-1][-1]
        schedule += schedule_c

    return schedule


def draw_polys(ps, ax):
    for p in ps:
        p_np = np.asarray(p.verts)
        poly = patches.Polygon(p_np, fill=False, edgecolor='0.5')
        ax.add_patch(poly)

    for p in ps:
        p_np = np.asarray(p.verts)
        plt.scatter(p_np[::-1, 0], p_np[::-1, 1])

    #ax.set_ylim(bottom=-1, top=6)
    #plt.show()
    plt.draw()



def get_curr_point(verts, t, speed):
    total_time = 0
    for i in range(len(verts)-1):
        p1 = np.array(verts[i]).astype(float)
        p2 = np.array(verts[i+1]).astype(float)
        total_time_ = total_time + l2_dist(verts[i], verts[i+1]) / speed
        if total_time_ > t:
            ratio = (t - total_time) / (total_time_ - total_time)
            new_p = (ratio*(p2 - p1)) + p1
            return new_p

        total_time = total_time_

    return verts[-1]


def animate_schedule(polygons, chedule, speed):
    fig, ax = plt.subplots()
    draw_polys(polygons, ax)
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

            if list(a) == s[-1][0] and list(b) == s[-1][1]:
                break
            else:
                #plt.clf()
                abc=1


if __name__ == '__main__':
    plt.ion()
    abc = 1
    speed = 15
    polygon = [[0, 0], [10, 0], [10, 10], [20, 20], [17, 28], [15, 30], [7, 15], [3, 15], [-5, 30], [-10, 20], [0, 10]]
    agents = [[0, 0], [10, 0]]

    # polygon = [[0, 0], [37, 0], [37, 8], [36, 11], [36, 22], [33, 22], [28, 24], [27, 20], [31, 21],
    #           [30, 19], [20, 19], [22, 24], [14, 16], [20, 14], [30, 14],
    #           [32, 11], [0, 10]]
    # agents = [[0, 0], [0, 10]]

    #polygon = [[0, 0], [40, 0], [40, 10], [38, 10], [38, 20]]
    #agents = [[0, 10], [0, 0]]
    # polygon = [[40, 0], [0, 60], [10, 50], [22, 55], [20, 45], [30, 42], [27, 50],
    #            [40, 40], [53, 50], [50, 42], [60, 45], [58, 55], [70, 50], [80, 60]]
    # polygon = polygon[::-1]
    # agents = [[10, 50], [0, 60]]

    # polygon = [[0,0], []]
    # polygon = polygon[::-1]
    # agents = [[10, 50], [0, 60]]

    #fig, ax = plt.subplots()
    #draw_polys([PolyNode(polygon)], ax)

    polygons = poly_decomp_cgal(polygon)
    vis_graph = create_vis_graph(polygon)

    dijkstra = DijkstraSPF(vis_graph, 1)
    path = dijkstra.get_path(4)

    polygons = [PolyNode(poly) for poly in polygons]
    #draw_polys(polygons, ax)

    create_edge(polygons)
    root = find_root_node(polygons, agents)
    schedule = dfs(polygons, root, polygon, vis_graph, [], agents)
    animate_schedule(polygons, schedule, speed)
    plt.pause(0)
    abc = 1