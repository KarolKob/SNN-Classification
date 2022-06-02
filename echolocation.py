from numpy import square
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

# Returns the closest point from the list to the point in the 2nd parameter
def closest_point(list, point):
    distance = 9999
    index = 0
    for i in range(0, len(list)):
        if list[i].distance(point) < distance:
            distance = list[i].distance(point)
            index = i
            
    return list[index]

# Find the side of a polygon that contains the intersection point and reflect the given line
def find_reflection(sides, l, inter):
    for ind in range(0, len(sides)):
        side = sides[ind]
        if (side.points[0][1] == inter[1] and side.points[1][1] == inter[1]) or (side.points[0][0] == inter[0] and side.points[1][0] == inter[0]):
            perpendicular = Line2D(side.points[0], side.points[1]).perpendicular_line(inter)
            print(side)
            print(perpendicular)
            sym_line = l.reflect(perpendicular)
            return sym_line
    return false

# Return the distance from the start of generating the sound to reflecting back
def count_reflections(circle_points, inner_polygon, outer_polygon):
    point_matrix = []
    for ind in range(0, len(circle_points)):
        circ = Circle(circle_points[ind], 0.25)
        rotating_point = Point2D(0, 0)
        pos_inc = 1
        neg_inc = -1

        # Repeat for different angles (ac)
        for rot_num in range(20):
            if rot_num % 2 == 1 and rot_num != 0:
                new_point = rotating_point.rotate((pi/180)*pos_inc, circle_points[ind])     # increasing 1 degree - pi/180
                pos_inc += 1
            elif rot_num % 2 == 0 and rot_num != 0:
                new_point = rotating_point.rotate((pi/180)*neg_inc, circle_points[ind])     # decreasing 1 degree
                neg_inc -= 1
            else:
                new_point = rotating_point

            l = Line2D(circle_points[ind], new_point)
            point_array = []
            point_array.append(circle_points[ind])

            inter = closest_point(inner_polygon.intersection(l), circle_points[ind])
            point_array.append(inter)

            # Find the side that contains the intersection point and reflect the line
            sym_line = find_reflection(inner_polygon.sides, l, inter)
            
            if len(circ.intersection(sym_line)) == 0:                                   # TODO: the intersection function can't count with larger complexity
                # Find the correct point of intersection with the big square
                inter = closest_point(outer_polygon.intersection(sym_line), inter)
                point_array.append(inter)

                sym_line2 = find_reflection(outer_polygon.sides, sym_line, inter)
                print("sym_line2: $", sym_line2)
                # Iterate until limit reached or intersected with the sensor
                for i in range(0, 10):
                    inter = closest_point(inner_polygon.intersection(sym_line2), inter)
                    point_array.append(inter)
                    sym_line = find_reflection(inner_polygon.sides, sym_line2, inter)
                    

                    if len(circ.intersection(sym_line)) > 0:
                        point_array.append(circle_points[ind])
                        point_matrix.append(point_array)
                        break
                    else:
                        inter = closest_point(outer_polygon.intersection(sym_line), inter)
                        point_array.append(inter)

                        sym_line2 = find_reflection(outer_polygon.sides, sym_line, inter)
            else:
                point_array.append(circle_points[ind])
                point_matrix.append(point_array)

# Virtual sources method for simulation of indoor acoustics

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

# Count the 2nd half of the circle
while x > -radius:
    y_square = pow(radius, 2) - pow(x, 2)

    if y_square == 0:
        y = 0
    else:
        y = -sqrt(y_square)

    cir_point = Point2D(x, y)
    cir_point_array.append(cir_point)
    x = x - 0.4
    count = count + 1

# Iterate through the circle creating segments from it to (0, 0)
l = Line2D(Point(0, 0), Point(2, 5))

inter = closest_point(sq.intersection(l), Point(2, 5))
print(inter)

sides = sq.sides

# Find the side that contains the intersection point and reflect the line
sym_line = find_reflection(sides, l, inter)

print(sym_line)

# Find the correct point of intersection with the big square
interb = closest_point(sqb.intersection(sym_line), inter)

print(interb)

sym_line2 = find_reflection(sqb.sides, sym_line, interb)

square_reflections = count_reflections(cir_point_array, sq, sqb)
print("Square reflections: $", square_reflections)

# l.reflect(Line()) reflect symmetrically
# Rotating a point against another point
#z.rotate(pi/4, inter)

#ang = l.angle_between(l2)
#print(ang)
