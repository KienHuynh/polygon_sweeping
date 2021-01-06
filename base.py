"""
This file provides basic functions to check for line segment intersection
and for checking the validity of a diagonal in a simple polygon
"""

def area2(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])


def left(a, b, c):
    return area2(a, b, c) > 0


def left_on(a, b, c):
    return area2(a, b, c) >= 0


def colinear(a, b, c):
    if a == b or a == c or b == c:
        return True
    return area2(a, b, c) == 0


# def inCone(polygon, a, b):
#     n = len(polygon)
#     a0 = polygon[(polygon.index(a) - 1 + n) % n]
#     a1 = polygon[(polygon.index(a) + 1) % n]
#
#     if left_on(a, a1, a0):
#         return left(a, b, a0) and left(b, a, a1)
#
#     return not (left_on(a, b, a1) and left_on(b, a, a0))


def xor(a, b):
    return (a and not b) or (not a and b)


def between(a, b, c):
    if not colinear(a, b, c):
        return False;

    if a[0] != b[0]:
        return ((a[0] <= c[0]) and (c[0] <= b[0])) or ((a[0] >= c[0]) and (c[0] >= b[0]))
    else:
        return ((a[1] <= c[1]) and (c[1] <= b[1])) or ((a[1] >= c[1]) and (c[1] >= b[1]))


def intersectProp(a, b, c, d):
    if (colinear(a, b, c) or
            colinear(a, b, d) or
            colinear(c, d, a) or
            colinear(c, d, b)):
        return False

    return xor(left(a, b, c), left(a, b, d)) and xor(left(c, d, a), left(c, d, b))


def intersect(a, b, c, d):
    if intersectProp(a, b, c, d):
        return True
    elif between(a, b, c) or between(a, b, d) or between(c, d, a) or between(c, d, b):
        return True
    return False


def diagonalie(a, b, verts):
    c = verts[0]
    i = 0
    a = verts[a]
    b = verts[b]
    for i in range(len(verts)):
        c = verts[i]
        c1 = verts[(i+1) % len(verts)]
        if ((c != a) and (c1 != a) and (c != b) and (c1 != b) and intersect(a, b, c, c1)):
            return False

    return True


def in_cone(a, b, verts):
    a1 = verts[(a + 1) % len(verts)]
    a0 = verts[(a - 1) % len(verts)]
    a = verts[a]
    b = verts[b]

    if (left_on(a, a1, a0)):
        return left(a, b, a0) and left(b, a, a1)

    return not (left_on(a, b, a1) and left_on(b, a, a0))


def diagonal(a, b, verts):
    return in_cone(a, b, verts) and in_cone(b, a, verts) and diagonalie(a, b, verts)
