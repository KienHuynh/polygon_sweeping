import matplotlib.pyplot as plt

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