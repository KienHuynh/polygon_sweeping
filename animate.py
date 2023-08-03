import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

from polygraph import *
from dijkstra import Graph, DijkstraSPF
from base import point_equal, l2_dist
import time


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


def draw_polygon(polygon, color='k'):
    n = len(polygon)
    polygon = np.asarray(polygon)
    plt.plot(polygon[:, 0], polygon[:, 1], color=color)
    plt.plot(polygon[[n-1, 0], 0], polygon[[n-1, 0], 1], color=color)


def draw_filled_polygon(polygon, color=[1, 0.9, 0.6], hatch=''):
    n = len(polygon)
    polygon = np.asarray(polygon)
    polygon = np.vstack((polygon, polygon[-1, :].reshape(1,2)))
    plt.fill(polygon[:, 0], polygon[:, 1], color = color, hatch=hatch)


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


def draw_poly_nodes(ps, ax, linewidth=0.5, edgecolor='0.7'):
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


def animate_schedule(polygons, subpolygons, schedule, speed, fig=None, ax=None, subpoly_color=[0.6, 0.6, 1]):
    """
    Animate the computed sweeping strategy
    :param polygon: List[List[float]], the original polygon and its holes
    :param polygons: List[PolyNode], list of convex subpolygons
    :param schedule: List[List[List[List[float]]]], the strategy.
    :param speed: float, animation speed
    :return: None
    """
    if not fig:
        fig, ax = plt.subplots()

    if type(subpolygons[0]) == PolyNode:
        draw_poly_nodes([PolyNode(p) for p in polygons], ax, edgecolor=[0, 0, 0], linewidth=1)
        draw_poly_nodes(subpolygons, ax)
    else:
        for sp in subpolygons:
            draw_polygon(sp, subpoly_color)

        for p in polygons:
            draw_polygon(p, 'k')

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


def animate_schedule_gif(polygons, subpolygons, schedule, speed, boundary_polygon, fig=None, ax = None, subpoly_color=[0.6, 0.6, 1]):
    """
    Animate the computed sweeping strategy
    :param polygon: List[List[float]], the original polygon and its holes
    :param polygons: List[PolyNode], list of convex subpolygons
    :param schedule: List[List[List[List[float]]]], the strategy.
    :param speed: float, animation speed
    :return: None
    """
    if fig == None:
        fig, ax = plt.subplots()
    for sp in subpolygons:
        draw_polygon(sp, subpoly_color)
    plt.plot(boundary_polygon[:, 0], boundary_polygon[:, 1], 'k')
    n = len(boundary_polygon)
    plt.plot(boundary_polygon[[n- 1, 0], 0], boundary_polygon[[n- 1, 0], 1], 'k')
    # rectangle = np.min(boundary_polygon, 0)
    # r_width, r_height = np.max(boundary_polygon, 0) - np.min(boundary_polygon, 0)
    # patch = patches.Rectangle(rectangle, r_width, r_height, facecolor = [1,1,1])
    # ax.add_patch(patch)
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
    segment_path = np.concatenate((segment_path, segment_path), 0)
    for s in schedule:
        t0 = time.time()
        path_a = [x[0] for x in s]
        path_b = [x[1] for x in s]
        xx = None
        while True:
            t1 = time.time()
            ela_time = t1 - t0
            a = get_curr_point(path_a, ela_time, speed)
            b = get_curr_point(path_b, ela_time, speed)
            segment_path[1, :, :] = np.asarray([[a[0], b[0]], [a[1], b[1]]])

            x = [segment_path[0, 0, 0], segment_path[0, 0, 1], segment_path[1, 0, 1], segment_path[1, 0, 0]]
            y = [segment_path[0, 1, 0], segment_path[0, 1, 1], segment_path[1, 1, 1], segment_path[1, 1, 0]]
            xx = plt.fill(x, y, facecolor = [1, 0.8, 0.5], alpha=1, edgecolor=[1, 0.8, 0.5])
            line.set_data([a[0], b[0]], [a[1], b[1]])
            scatter.set_offsets([[a[0], a[1]], [b[0], b[1]]])
            plt.show()
            plt.axis('equal')
            plt.pause(0.00000001)

            segment_path[0, :, :] = segment_path[1, :, :]

            try:
                if point_equal(a, s[-1][0]) and point_equal(b, s[-1][1]): #list(a) == s[-1][0] and list(b) == s[-1][1]:
                    break
            except Exception as e:
                abc = 1