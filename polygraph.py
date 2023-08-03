from typing import List
import numpy as np

from animate import *
from base import *
import copy

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



class SubPolygon:
    def __init__(self, polygon: List[List], polytype: str, base_edge=None, left_edge=None, right_edge=None):
        """
        :param polygon: [[x1, y1], [x2, y2], ...]
        :param type: "hist" or "vis"
        :param base_edge: [[x1, y1], [x2, y2]] base edge of the histogram
        :param left_edge: [[x1, y1], [x2, y2]] left edge of the histogram (oriented upward), can be a single point
        :param right_edge: [[x1, y1], [x2, y2]] right edge of the histogram (oriented upward), can be a single point
        """
        self.polygon = polygon
        self.type = polytype
        if self.type == 'hist':
            self.base_edge =  base_edge

        self.right_edge = right_edge
        self.left_edge = left_edge


    def compute_hist_schedule(self, starting_point=False, ending_point=False):
        """

        :param starting_point: if True, assume that the starting point is the bottom right vertex of the histogram,
        then add the transition from this starting point to the right edge
        :return: schedule, t
        """
        origin = self.base_edge[0]
        unit_v = np.asarray([1, 0])
        base_edge = np.asarray(self.base_edge) - origin
        polygon = np.asarray(self.polygon)
        polygon = polygon[2:, :]
        projections = np.copy(polygon) - origin
        cos_val = get_cos(base_edge[1] - base_edge[0], unit_v)
        # Get the angle from Ox to base_edge vector
        if base_edge[1, 1] < 0:
            sin_val = 0 - np.sqrt(1 - cos_val**2)
        else:
            sin_val = np.sqrt(1 - cos_val ** 2)

        # Project all points of the upper chain of the histogram polygon onto the base edge
        r_mat_ox = np.asarray([[cos_val, sin_val], [-sin_val, cos_val]])
        r_mat = np.asarray([[cos_val, -sin_val], [sin_val, cos_val]])
        projections = (r_mat_ox@projections.T).T
        projections[:, 1] = 0
        projections = (r_mat@projections.T).T + origin

        # Create the actual schedule
        #polygon = polygon[::-1, :]
        #projections = projections[::-1, :]
        schedule = []
        if starting_point:
            schedule += [
                [
                    [self.right_edge[0], self.right_edge[0]],
                    self.right_edge
                ]
            ]

        schedule += [
                    [
                        self.right_edge,
                        [list(projections[0, :]), list((polygon[0, :]))]
                    ]
                ]

        n = len(polygon)
        for i in range(n-1):
            s = [
                [list(projections[i, :]), list(polygon[i, :])],
                [list(projections[(i + 1), :]), list(polygon[(i + 1), :])]
                ]
            schedule.append(s)

        s = [
            [list(projections[n-1, :]), list(polygon[n-1, :])],
            list(self.left_edge)
        ]
        schedule.append(s)

        if ending_point:
            schedule += [
                [
                    self.left_edge,
                    [self.left_edge[1], self.left_edge[1]]
                ]
            ]

        return schedule


    def compute_vis_schedule(self, starting_point=False, ending_point=False):
        schedule = []
        origin = self.polygon[0]
        if starting_point:
            schedule = [
                [
                    [self.right_edge[0], self.right_edge[0]],
                    self.right_edge
                ]
            ]
        for i in range(1, len(self.polygon)-1):
            s = [
                [origin, list(self.polygon[i])],
                [origin, list(self.polygon[(i + 1)])]
                ]
            schedule.append(s)

        if ending_point:
            schedule += [
                [
                    self.left_edge,
                    [self.left_edge[0], self.left_edge[0]]
                ]
            ]
        return schedule


    def compute_schedule(self, starting_point=False, ending_point=False):
        """
        If histogram:
        assume the starting configuration is the right edge, sweep until the segment becomes the left edge
        If vis:
        assume the starting configuration is the right edge, (radially) sweep until the segment becomes the left edge
        :return:
        """
        schedule = []
        if self.type == 'hist':
            schedule = self.compute_hist_schedule(starting_point, ending_point)
        if self.type == 'vis':
            schedule = self.compute_vis_schedule(starting_point, ending_point)

        t = 0
        for s in schedule:
            for i in range(len(s) - 1):
                e0 = s[i]
                e1 = s[i + 1]
                da = l2_dist(e0[0], e1[0])
                db = l2_dist(e0[1], e1[1])
                t += max(da, db)
        return schedule, t


class HistogramNode:
    def __init__(self):
        # This stores the actual polygons from the vis/histogram decomposition
        # Note that it is ordered ('right' to 'left')
        self.subpolygons = []
        self.dividers = []  # These are the chords inbetween the subpolygons

        # Children of this node, ordered clockwise around the boundary of this node
        self.children = []

        # The outer boundary of this node (boundary of the union of the subpolygons)
        self.outer_boundary = []

        # Visibility graph of the polygon of the outer boundary
        self.vis_graph = []

        # This is the interface edge
        self.interface_edge = []

    def compute_outer_boundary(self):

        # Get all of the base edges of the histogram, i.e. the bottom reflex chain
        last_p = []
        for p in self.subpolygons[::-1]:
            if p.type == 'hist':
                self.outer_boundary.append(p.base_edge[0])
                last_p = p

        # Get the top polygonal chain
        self.outer_boundary.append(last_p.base_edge[1])

        last_p = []
        for p in self.subpolygons:
            bdy = []
            if p.type == 'vis':
                bdy += p.polygon[1:-1]
            if p.type == 'hist':
                bdy += p.polygon[2:-1]
            self.outer_boundary += bdy
            last_p = p
        self.outer_boundary.append(last_p.polygon[-1])

        for i in range(len(self.outer_boundary)):
            j = (i + 1) % len(self.outer_boundary)
            pi = np.asarray(self.outer_boundary[i])
            pj = np.asarray(self.outer_boundary[j])
            #assert np.sum(pi - pj) != 0, "Outer boundary has repeating vertices"

        self.vis_graph = create_vis_graph(self.outer_boundary, [])
        # self.outer_boundary = list(np.unique(self.outer_boundary, axis=0))


    def compute_shortest_path(self, s, t):
        """
        Using visibility graph
        :param s:
        :param t:
        :return:
        """
        si = [i for i in range(len(self.outer_boundary)) if point_equal(self.outer_boundary[i], s)][0]
        ti = [i for i in range(len(self.outer_boundary)) if point_equal(self.outer_boundary[i], t)][0]
        dijkstra = DijkstraSPF(self.vis_graph, si)
        path = dijkstra.get_path(ti)
        return path


    def get_all_outer_boundary(self):
        result = [self.outer_boundary]
        for c in self.children:
            result += c.get_all_outer_boundary()
        return result

    def compute_schedule(self, no_return=False, next_sibling_point=None, visualize=False):
        schedule = []
        time = 0

        # Compute the sweeping schedule for the histogram polygon and the vis polygons of this node
        for i in range(len(self.subpolygons)):
            p = self.subpolygons[i]
            starting_point = False
            if i == 0:
                starting_point = True
            s, t = p.compute_schedule(starting_point=starting_point)
            schedule += s
            time += t
        current_edge = schedule[-1][-1]

        for i in range(len(self.children)):
            # Move to child node
            c = self.children[i]
            ie = c.interface_edge
            stp_a = self.compute_shortest_path(current_edge[0], ie[1])
            stp_b = self.compute_shortest_path(current_edge[1], ie[1])
            stp_a = [self.outer_boundary[x] for x in stp_a]
            stp_b = [self.outer_boundary[x] for x in stp_b]
            if len(stp_b) == len(stp_a)-1:
                stp_b.insert(0, stp_b[0])
            elif len(stp_a) == len(stp_b)-1:
                stp_a.insert(0, stp_a[0])
            s = [[current_edge, [stp_a[0], stp_b[0]]]]
            for j in range(len(stp_a)-1):
                s.append([
                    [stp_a[j], stp_b[j]],
                    [stp_a[j+1], stp_b[j+1]]
                ])

            schedule += s

            # Compute the schedule of the child node
            if i == len(self.children)-1:
                c_no_return = True
            else:
                c_no_return = False
            c_schedule = c.compute_schedule(no_return and c_no_return)
            schedule += c_schedule
            current_edge = c_schedule[-1][-1]

        if no_return:
            return schedule

        ie = self.interface_edge

        stp_a = self.compute_shortest_path(current_edge[0], ie[1])
        stp_b = self.compute_shortest_path(current_edge[1], ie[1])
        stp_a = [self.outer_boundary[x] for x in stp_a]
        stp_b = [self.outer_boundary[x] for x in stp_b]
        if len(stp_b) == len(stp_a) - 1:
            stp_b.insert(0, stp_b[0])
        elif len(stp_a) == len(stp_b) - 1:
            stp_a.insert(0, stp_a[0])

        s = [[current_edge, [stp_a[0], stp_b[0]]]]
        for i in range(len(stp_a) - 1):
            s.append([
                [stp_a[i], stp_b[i]],
                [stp_a[i + 1], stp_b[i + 1]]
            ])
        schedule += s

        return schedule