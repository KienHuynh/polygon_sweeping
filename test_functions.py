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
        [208, 680], [200, 680], [200, 688], [196, 692], [192, 684]]

    polygon1 = np.asarray(polygon1, dtype=np.float64)
    n1 = len(polygon1)
    mat = np.asarray([[-0.935458, -0.353438], [0.353438, -0.935458]]).T
    polygon1 = (mat@polygon1.T).T
    base_edge1 = polygon1[[0, 1], :]
    left_edge1 = polygon1[[0, n1-1], :]
    right_edge1= polygon1[[1, 2], :]
    speed1 = 40

    polygon = polygon1
    base_edge = base_edge1
    left_edge = left_edge1
    right_edge = right_edge1
    speed = speed1

    t = 0
    for i in range(2, len(polygon) - 1):
        a = polygon[i]
        b = polygon[i + 1]
        t += l2_dist(a, b)
    print(t)

    hist = SubPolygon(polygon, 'hist', base_edge, left_edge, right_edge)
    schedule, t = hist.compute_schedule()
    print(t)
    plt.ion()
    # fig, ax = plt.subplots()
    # for s in schedule[0]:
    #
    #     draw_polygon(polygon)
    #     plt.plot([s[0][0], s[1][0]], [s[0][1], s[1][1]], color=[0.7, 0.7, 0.7])
    #     plt.show()
    #     plt.axis('equal')
    #     plt.waitforbuttonpress()
    # plt.close(fig)

    animate_schedule([list(polygon)], [], schedule, speed)
    plt.waitforbuttonpress()
    plt.close()


def test_vis_sweep():
    polygon = [[216.134, 601.741], [219.497, 584.058], [235.228, 588.18], [225.247, 600.548], [238.699, 607.491],
        [228.285, 608.576], [238.265, 618.665], [236.638, 627.887], [230.237, 632.986], [223.945, 624.09],
        [216.242, 631.792], [197.04, 628.755], [206.262, 615.519], [190.314, 618.774], [204.526, 602.718],
        [189.989, 599.571], [193.135, 592.303], [208.106, 597.727], [208.323, 591.652]]
    left_edge = [[216.134, 601.741], [208.323, 591.652]]
    right_edge = [[216.134, 601.741], [219.497, 584.058]]
    speed = 20
    t = 0
    for i in range(1, len(polygon) - 1):
        a = polygon[i]
        b = polygon[i + 1]
        t += l2_dist(a, b)
    print(t)

    hist = SubPolygon(polygon, 'vis', left_edge=left_edge, right_edge=right_edge)
    schedule, t = hist.compute_schedule()
    print(t)
    plt.ion()
    animate_schedule([list(polygon)], [], schedule, speed)
    plt.waitforbuttonpress()
    plt.close()


if __name__ == '__main__':
    # test_vis()
    # test_hist()
    test_hist_sweep()
    test_vis_sweep()