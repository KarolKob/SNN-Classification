from sympy import *
from sympy.geometry import *
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import numpy as np

# Returns the closest point from the list to the point in the 2nd parameter
def closest_point(list, point):
    distance = 9999
    index = 0
    for i in range(0, len(list)):
        new_dist = list[i].distance(point)
        if new_dist < distance and new_dist != 0:
            distance = new_dist
            index = i
            
    return list[index]

def closest_point_same_side(list, point, side, prev_point):
    distance = 9999
    index = 0
    A = side.p1
    B = side.p2
    prev_position = sign((B.x - A.x) * (prev_point.y - A.y) - (B.y - A.y) * (prev_point.x - A.x))
    for i in range(0, len(list)):
        position = sign((B.x - A.x) * (list[i].y - A.y) - (B.y - A.y) * (list[i].x - A.x))

        if position == prev_position:
            new_dist = list[i].distance(point)
            if new_dist < distance and new_dist != 0:
                distance = new_dist
                index = i
    
    if distance < 200:
        return list[index]
    else:
        return False

# Find the side of a polygon that contains the intersection point and reflect the given line
def find_reflection(sides, l, inter):
    #print("inter: ", inter)
    #print("sides: ", sides)
    for ind in range(0, len(sides)):
        side = sides[ind]
        if side.distance(inter) == 0:
            perpendicular = Line2D(side.points[0], side.points[1]).perpendicular_line(inter)
            sym_line = l.reflect(perpendicular)
            #print("reflection found: ", sym_line)
            return sym_line, side
    return False

class IntersectionPlus:
    def __init__(self, middle_point, line_len) -> None:
        self.segHor = Segment2D(Point2D(middle_point.x - line_len, middle_point.y), 
                                    Point2D(middle_point.x + line_len, middle_point.y))
        self.segVer = Segment2D(Point2D(middle_point.x, middle_point.y - line_len), 
                                    Point2D(middle_point.x, middle_point.y + line_len))

    def intersection_found(self, line):
        inter = self.segHor.intersection(line)
        if inter != []:
            return True
        else:
            inter = self.segVer.intersection(line)
            if inter != []:
                return True
            else:
                return False

def pick_obj_and_find_reflection(big_obj, small_obj, to_reflect, inter, side, prev_inter):
    inter_list = small_obj.intersection(to_reflect)

    if inter_list != []:
        inter_s = closest_point_same_side(inter_list, inter, side, prev_inter)
        inter_b = closest_point_same_side(big_obj.intersection(to_reflect), inter, side, prev_inter)

        if  inter_s != False and inter_s.distance(inter) < inter_b.distance(inter):
            return find_reflection(small_obj.sides, to_reflect, inter_s), inter_s
        else:
            return find_reflection(big_obj.sides, to_reflect, inter_b), inter_b
    else:
        inter_b = closest_point_same_side(big_obj.intersection(to_reflect), inter, side, prev_inter)
        return find_reflection(big_obj.sides, to_reflect, inter_b), inter_b
        


# Return the distance from the start of generating the sound to reflecting back
def angle_loop(ind, circle_points, inner_polygon, outer_polygon):
    point_matrix = []
    all_index = 0

    rotating_point = Point2D(0, 0)
    pos_inc = 1
    neg_inc = -1

    inter_plus = IntersectionPlus(circle_points[ind], 0.1)

    # Repeat for different angles (ac)
    for rot_num in range(20):
        if rot_num % 2 == 1 and rot_num != 0:
            new_point = rotating_point.rotate((pi/180)*pos_inc, circle_points[ind])     # increasing 1 degree - pi/180
            new_point = Point2D(new_point.x.evalf(), new_point.y.evalf())
            pos_inc += 1
        elif rot_num % 2 == 0 and rot_num != 0:
            new_point = rotating_point.rotate((pi/180)*neg_inc, circle_points[ind])     # decreasing 1 degree
            new_point = Point2D(new_point.x.evalf(), new_point.y.evalf())
            neg_inc -= 1
        else:
            new_point = rotating_point

        #all_index += 1
        #print(all_index)
        l = Line2D(circle_points[ind], new_point)
        point_array = []
        point_array.append(circle_points[ind])

        inter = closest_point(inner_polygon.intersection(l), circle_points[ind])
        point_array.append(inter)

        # Find the side that contains the intersection point and reflect the line
        sym_line, side = find_reflection(inner_polygon.sides, l, inter)
        
        if not inter_plus.intersection_found(sym_line):
            # Iterate until limit reached or intersected with the sensor
            for i in range(0, 10):
                refl_list, inter = pick_obj_and_find_reflection(outer_polygon, inner_polygon, sym_line, inter, side, point_array[len(point_array)-2])
                sym_line = refl_list[0]
                side = refl_list[1]
                point_array.append(inter)
                
                if inter_plus.intersection_found(sym_line):
                    point_array.append(circle_points[ind])
                    point_matrix.append(point_array)
                    break

                if i == 9:
                    point_matrix.append(point_array)
        else:
            point_array.append(circle_points[ind])
            point_matrix.append(point_array)

    min_dist_arr = []
    for i in range(0, len(point_matrix)):
        min_dist = 99999
        point_dist = 0
        for j in range(1, len(point_matrix[i])):
            point_dist += point_matrix[i][j].distance(point_matrix[i][j - 1])
        if min_dist > point_dist:
            min_dist = point_dist.evalf()
        min_dist_arr.append(min_dist)

    min_dist = 99999
    for i in range(0, len(min_dist_arr)):
        if min_dist > min_dist_arr[i]:
            min_dist = min_dist_arr[i]

    return min_dist       # Count the min value for each iteration

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

#Reflections of a square to check the network
square_reflections = Parallel(n_jobs=16)(delayed(angle_loop)(ind, cir_point_array, sq, sqb) for ind in range(0, len(cir_point_array)))
print("Square reflections: $", square_reflections)

np.savetxt('Square.txt', square_reflections, fmt='%.6f')

# Points for the letter N
n1, n2, n3, n4, n5 = [(1, 1), (0.5, 1), (0.5, -0.33), (-0.5, 1), (-1, 1)]
n6, n7, n8, n9, n10 = [(-1, -1), (-0.5, -1), (-0.5, 0.33), (0.5, -1), (1, -1)]
letterN = Polygon(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)

# Points for the letter A
a1, a2, a3, a4 = [(0.2, 1), (-0.2, 1), (-1, -1), (-0.6, -1)]
l1 = Line2D(Point(-1, -1), Point(-0.2, 1))
l2 = l1.parallel_line(Point(-0.6, -1))
a5 = l2.intersection(Line2D(Point(1, -0.5), Point(-1, -0.5)))[0]
a5, a6, a7, a8 = [(a5.x, a5.y), (-1*a5.x, a5.y), (0.6, -1), (1, -1)]
letterA = Polygon(a1, a2, a3, a4, a5, a6, a7, a8)

# Points for the letter X
x1, x2, x3, x4, x5, x6 = [(1, 1), (0.5, 1), (0, 0.2), (-0.5, 1), (-1, 1), (-0.3, 0)]
x7, x8, x9, x10, x11, x12 = [(-1, -1), (-0.5, -1), (0, -0.2), (0.5, -1), (1, -1), (0.3, 0)]
letterX = Polygon(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)

letterA_reflections = Parallel(n_jobs=16)(delayed(angle_loop)(ind, cir_point_array, letterA, sqb) for ind in range(0, len(cir_point_array)))
print("Letter A reflections: $", letterA_reflections)

np.savetxt('Letter_A.txt', letterA_reflections, fmt='%.6f')

letterN_reflections = Parallel(n_jobs=16)(delayed(angle_loop)(ind, cir_point_array, letterN, sqb) for ind in range(0, len(cir_point_array)))
print("Letter N reflections: $", letterN_reflections)

np.savetxt('Letter_N.txt', letterN_reflections, fmt='%.6f')

letterX_reflections = Parallel(n_jobs=16)(delayed(angle_loop)(ind, cir_point_array, letterX, sqb) for ind in range(0, len(cir_point_array)))
print("Letter X reflections: $", letterX_reflections)

np.savetxt('Letter_X.txt', letterX_reflections, fmt='%.6f')

#X = np.array([1, 0.5, 0.5, -0.5, -1, -1, -0.5, -0.5, 0.5, 1])
#Y = np.array([1, 1, -0.33, 1, 1, -1, -1, 0.33, -1, -1])

# Plotting point using scatter method
#plt.scatter(X,Y)
#plt.show()