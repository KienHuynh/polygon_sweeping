from two_sweeping_approx import *


def draw_polygon(polygon, color='k'):
    n = len(polygon)
    plt.plot(polygon[:, 0], polygon[:, 1], color)
    plt.plot(polygon[[n-1, 0], 0], polygon[[n-1, 0], 1], color)


def draw_filled_polygon(polygon, color=[1, 0.9, 0.6]):
    n = len(polygon)
    polygon = np.vstack((polygon, polygon[-1, :].reshape(1,2)))
    plt.fill(polygon[:, 0], polygon[:, 1], color = color)

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

    polygon6 = [[208, 560], [200, 552], [192, 560], [200, 564], [192, 572],
        [212, 572], [200, 584], [220, 580], [228, 584], [236, 576],
        [220, 568], [236, 556], [224, 560]]
    polygon6 = polygon6[::-1]
    polygon6 = polygon6[9:] + polygon6[0:9]

    polygon7 = [[224, 608], [240, 608], [240, 624], [256, 624], [256, 656], [240, 656],
        [240, 720], [224, 720], [224, 704], [216, 704], [216, 692], [224, 692],
        [224, 680], [216, 680], [216, 668], [224, 668], [224, 656], [204, 656],
        [204, 680], [184, 680], [184, 668], [196, 668], [196, 640], [176, 640],
        [176, 624], [224, 624]]

    polygon = polygon7
    for i in range(0, len(polygon)):
        edge_i = [polygon[i], polygon[(i + 1) % len(polygon)]]
        compute_histogram(polygon, edge_i)


def test_hist_sweep():
    polygon1 = [[192, 664], [260, 664], [260, 696], [248, 688], [248, 700], [240, 688],
        [228, 692], [228, 700], [224, 704], [224, 716], [216, 720], [208, 712],
        [208, 680], [200, 680], [200, 688], [196, 692], [192, 684], [192, 664]]
    polygon1 = np.asarray(polygon1, dtype=np.float64)
    n1 = len(polygon1)
    mat = np.asarray([[-0.709516, 0.704689], [-0.704689, -0.709516]]).T
    mat = np.asarray([[-0.709516, 0.704689], [-0.704689, -0.709516]]).T
    polygon1 = (mat@polygon1.T).T
    base_edge1 = polygon1[[0, 1], :]
    left_edge1 = polygon1[[0, n1-1], :]
    right_edge1= polygon1[[1, 2], :]


    polygon = polygon1
    base_edge = base_edge1
    left_edge = left_edge1
    right_edge = right_edge1
    hist = SubPolygon(polygon, 'hist', base_edge, left_edge, right_edge)
    upper_chain, projected_points = hist.computer_hist_schedule()
    fig, ax = plt.subplots()
    draw_polygon(polygon1)
    plt.scatter(projected_points[:, 0], projected_points[:, 1], c='r')
    for i in range(len(upper_chain)):
        plt.plot([upper_chain[i, 0], projected_points[i, 0]], [upper_chain[i, 1], projected_points[i, 1]], c=[0.7, 0.7, 0.7])

    plt.axis('equal')
    plt.show()
    plt.close(fig)


if __name__ == '__main__':
    # test_vis()
    # test_hist()
    test_hist_sweep()