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

pb1, pb2, pb3, pb4 = [(0, 4), (-4, 0), (0, -4), (4, 0)]
sqb = Polygon(pb1, pb2, pb3, pb4)
print(sqb.area)


radius = 2
x = -2
cir_point_array = []
sampling_param = 0.4    # 0.4 for 20 samples overall

print("Circle time")
count = 0
# Circle function y^2 = r^2 - x^2
# Count the first half of the sensor trajectory
while x < radius:
    y_square = pow(radius, 2) - pow(x, 2)

    if y_square == 0:
        y = 0
    else:
        y = sqrt(y_square)

    cir_point = Point(x, y)
    cir_point_array.append(cir_point)
    x = x + 0.4
    count = count + 1

print(count)

x = radius

while x > -radius:
    y_square = pow(radius, 2) - pow(x, 2)

    if y_square == 0:
        y = 0
    else:
        y = -sqrt(y_square)

    cir_point = Point(x, y)
    cir_point_array.append(cir_point)
    x = x - 0.4
    count = count + 1

print(count)
print(cir_point_array)

# Virtual sources method for simulation of indoor acoustics