from sympy import *
from sympy.geometry import *
x = Point(0, 0)
y = Point(1, 1)
z = Point(2, 2)
zp = Point(1, 0)
t = Triangle(zp, y, x)
print(t.area)
print(t.medians[x])
m = t.medians
print(intersection(m[x], m[y], m[zp]))
c = Circle(x, 5)
l = Line(Point(5, -5), Point(5, 5))
print(c.is_tangent(l)) # is l tangent to c?
l = Line(x, y)
print(c.is_tangent(l)) # is l tangent to c?
print(intersection(c, l))
p1, p2, p3, p4 = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
sq = Polygon(p1, p2, p3, p4)
print(sq.area)

# Virtual sources method for simulation of indoor acoustics