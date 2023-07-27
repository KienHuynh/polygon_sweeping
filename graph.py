from typing import List
import numpy as np

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
        self.polygon = np.asarray(polygon, dtype=np.float64)
        self.type = polytype
        if self.type == 'hist':
            self.base_edge = np.asarray(base_edge, dtype=np.float64)

        self.right_edge = right_edge
        self.left_edge = left_edge


    def compute_hist_schedule(self):
        origin = self.base_edge[0]
        unit_v = np.asarray([1, 0])
        base_edge = self.base_edge - origin
        polygon = np.copy(self.polygon)
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
        schedule = [
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
        return schedule


    def compute_vis_schedule(self):
        schedule = []
        origin = list(self.polygon[0, :])
        for i in range(1, len(self.polygon)-1):
            s = [
                [origin, list(self.polygon[i, :])],
                [origin, list(self.polygon[(i + 1), :])]
                ]
            schedule.append(s)
        return schedule


    def compute_schedule(self):
        """
        If histogram:
        assume the starting configuration is the right edge, sweep until the segment becomes the left edge
        If vis:
        assume the starting configuration is the right edge, (radially) sweep until the segment becomes the left edge
        :return:
        """
        schedule = []
        if self.type == 'hist':
            schedule = self.compute_hist_schedule()
        if self.type == 'vis':
            schedule = self.compute_vis_schedule()

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
        # Note that it is ordered
        self.subpolygons = []
        self.dividers = []  # These are the chords inbetween the subpolygons

        # Children of this node, ordered clockwise around the boundary of this node
        self.children = []

        # The outer boundary of this node (boundary of the union of the subpolygons)
        self.outer_boundary = []

        # This is the interface edge[]