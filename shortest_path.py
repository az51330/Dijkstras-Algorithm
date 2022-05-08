"""
    File:        shortest_path.py
    Author:      A.J. Zuckerman
    Course:      CS 307 - Computational Geometry
    Assignment:  Synthesis - Shortest Path in Polygons
    Description: Finds the shortest path from start to goal point
    in a polygon with obstacles using Dijkstra's algorithm, runs experiment of
    query time to show O(nlogn + k) time, builds structure brute force since not testing that
"""

import math, random, heapq
import matplotlib.pyplot as plt
from matplotlib import collections as mpl

class Vertex:
    def __init__(self, x, y):
        self.loc = (x, y)
        self.cost = 0
        # an arc is a connection of two visible vertices listed by self loc then connect loc
        self.arcs = []
        self.weights = {}
        self.path = []

    def __str__(self):
        return f"({self.loc[0]}, {self.loc[1]})"

    def __lt__(self, other):
        return self.cost < other.cost

def euclidean_distance(point1, point2):
    """
    Calculates the euclidean distance of two given points
    
    keyword arguments:
    point1 -- tuple of x and y
    point2 -- tuple of x and y

    returns the float distance of the two points
    """
    return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)

def orient(p, q, r):
    """
    Orient calculates the orientation of three points using determinant method

    keyword arguments: 
    p -- start point of the edge
    q -- end point of the edge
    r -- point in question

    returns 1 if left, -1 if right, and 0 if straight
    """
    return q[0] * r[1] + p[0] * q[1] + r[0] * p[1] - q[0] * p[1] - r[0] * q[1] - p[0] * r[1]

def generate_quad(n):
    """
    Helper for generate obstacles, generates n^2 obstacles mapping

    keyword arguments:
    n -- number of vertices

    returns a list of points and edges
    """
    bottom = n // 2
    mid = bottom // 2
    mid_bottom = Vertex(mid, -(1 / bottom) * mid * (mid - bottom))
    mid_top = Vertex(mid, (1 / bottom) * mid * (mid - bottom) + bottom)
    points = [mid_top, mid_bottom]
    edges = []
    last_bottom_right = last_bottom_left = points[0]
    last_top_right = last_top_left = points[1]

    for index in range(mid):
        new_bottom_right = Vertex(mid + index, -(1 / bottom) * (mid + index) * (mid + index - bottom))
        points.append(new_bottom_right)
        edges.append((last_bottom_right, new_bottom_right))
        last_bottom_right = new_bottom_right

        new_bottom_left = Vertex(mid - index, -(1 / bottom) * (mid - index) * (mid - index - bottom))
        points.append(new_bottom_left)
        edges.append((last_bottom_left, new_bottom_left))
        last_bottom_left = new_bottom_left

        new_top_right = Vertex(mid + index, (1 / bottom) * (mid + index) * (mid + index - bottom) + bottom)
        points.append(new_top_right)
        edges.append((last_top_right, new_top_right))
        last_top_right = new_top_right

        new_top_left = Vertex(mid - index, (1 / bottom) * (mid - index) * (mid - index - bottom) + bottom)
        points.append(new_top_left)
        edges.append((last_top_left, new_top_left))
        last_top_left = new_top_left

    edges.append((last_bottom_left, last_bottom_right))
    edges.append((last_top_left, last_top_right))
    edges.append((last_bottom_left, mid_bottom))
    edges.append((last_bottom_right, mid_bottom))
    edges.append((last_top_left, mid_top))
    edges.append((last_top_right, mid_top))

    edges.append((Vertex(mid - .1, -(1 / bottom) * (mid - .1) * (mid - .1 - bottom)), Vertex(mid - .1, 1)))
    edges.append((Vertex(mid + .1, -(1 / bottom) * (mid + .1) * (mid + .1 - bottom)), Vertex(mid + .1, 1)))
    edges.append((Vertex(mid - .1, (1 / bottom) * (mid - .1) * (mid - .1 - bottom) + bottom), Vertex(mid - .1, bottom)))
    edges.append((Vertex(mid + .1, (1 / bottom) * (mid + .1) * (mid + .1 - bottom) + bottom), Vertex(mid + .1, bottom)))
    
    edges.append((Vertex(bottom - 1.1, -(1 / bottom) * (bottom - 1.1) * (bottom - 1.1 - bottom)), Vertex(bottom - 1.1, 1)))
    edges.append((Vertex(1 + .1, -(1 / bottom) * (1 + .1) * (1 + .1 - bottom)), Vertex(1 + .1, 1)))
    edges.append((Vertex(bottom - 1.1, (1 / bottom) * (bottom - 1.1) * (bottom - 1.1 - bottom) + bottom), Vertex(bottom - 1.1, bottom)))
    edges.append((Vertex(1 + .1, (1 / bottom) * (1 + .1) * (1 + .1 - bottom) + bottom), Vertex(1 + .1, bottom)))

    return points, edges

def generate_nlogn(n):
    """
    Helper for generate obstacles, generates nlogn obstacles mapping

    keyword arguments:
    n -- number of vertices

    returns a list of points and edges
    """
    points = []
    edges = []
    min_x_safe = 0
    min_y_safe = 0
    while n > 0:
        first = last = Vertex(random.randint(min_x_safe, min_x_safe + 5), random.randint(min_y_safe, min_y_safe + 5))
        points.append(first)
        for index in range(2):
            new_vert = Vertex(random.randint(min_x_safe, min_x_safe + 5), random.randint(min_y_safe, min_y_safe + 5))
            points.append(new_vert)
            edges.append((last, new_vert))
            last = new_vert
        edges.append((first, last))
        n -= 3
        min_x_safe += 6
        min_y_safe += 6

    return points, edges



def generate_obstacles(n, type):
    """
    Generates n^2 or nlogn obstacles mapping

    keyword arguments:
    n -- number of vertices
    type -- time complexity wanted

    returns a list of points and edges
    """
    if type == "quadratic":
        return generate_quad(n)
    elif type == "nlogn":
        return generate_nlogn(n)

def naive_visibility_graph(points, edges, start):
    """
    Creates a visibility graph by brute force in n^3 time

    keyword arguments:
    points -- a list of points in the visibility graph (points on obstacles)
    edges -- all edges of obstacles
    start -- the start vertex

    returns the start vertex which connects to the rest of the graph through its arcs
    """
    # CITE: http://www.science.smith.edu/~istreinu/Teaching/Courses/274/Spring98/Projects/Philip/fp/visibility.htm
    # Gave me pseudocode for the naive algorithm
    points = [start] + points
    for point in points:
        for other in points:
            arc = True
            if point != other:
                for edge in edges:
                    if edge == (point, other) or edge == (other, point):
                        arc = False
                    if orient(edge[0].loc, edge[1].loc, point.loc) < 0 and orient(edge[0].loc, edge[1].loc, other.loc) > 0 and orient(point.loc, other.loc, edge[0].loc) < 0 and orient(point.loc, other.loc, edge[1].loc) > 0:
                            arc = False
                    elif orient(edge[0].loc, edge[1].loc, point.loc) > 0 and orient(edge[0].loc, edge[1].loc, other.loc) < 0 and orient(point.loc, other.loc, edge[0].loc) < 0 and orient(point.loc, other.loc, edge[1].loc) > 0:
                            arc = False
                    elif orient(edge[0].loc, edge[1].loc, point.loc) < 0 and orient(edge[0].loc, edge[1].loc, other.loc) > 0 and orient(point.loc, other.loc, edge[0].loc) > 0 and orient(point.loc, other.loc, edge[1].loc) < 0:
                            arc = False
                    elif orient(edge[0].loc, edge[1].loc, point.loc) > 0 and orient(edge[0].loc, edge[1].loc, other.loc) < 0 and orient(point.loc, other.loc, edge[0].loc) > 0 and orient(point.loc, other.loc, edge[1].loc) < 0:
                            arc = False
                if arc:
                    point.arcs.append(other)
                    point.weights[other] = euclidean_distance(point.loc, other.loc)
    return start

def dijkstras_alg(points, edges, start, end):
    """
    computes dijkstra's algorithm on a visibility graph to determine the shortest path from a start to end point

    keyword arguments:
    points -- a list of points in the visibility graph (points on obstacles)
    edges -- all edges of obstacles
    start -- the start vertex
    end -- the end vertex

    returns a list of the arcs of the shortest path from the start to end vertices
    """
    queue = [start]
    heapq.heapify(queue)

    # Sets constant time look up hash tables might not need it though the algorithm should discriminate
    # against revisiting since path cost is added together
    explored = set()

    for point in points:
        if point != start:
            point.cost = math.inf

    while len(queue) > 0:
        curr = heapq.heappop(queue)
        explored.add(curr)
        if curr == end:
            final_path = curr.path + [curr]
            return final_path
        else:
            for arc in curr.arcs:
                if (curr.cost + curr.weights[arc]) < arc.cost and arc not in explored: 
                    heapq.heappush(queue, arc)
                    arc.path = curr.path + [curr]
                    arc.cost = curr.cost + curr.weights[arc]
    
    return []


def test_arcs(points, edges, path):
    """
    prints arcs of all points to check if the go through an obstacle

    keyword arguments:
    points -- a list of points in the visibility graph (points on obstacles)
    edges -- all edges of obstacles
    start -- the start vertex
    """
    # CITE: https://matplotlib.org/3.5.0/gallery/shapes_and_collections/line_collection.html
    # DESC: Never really used matplotlib too much so here is some line collection documentation that helped me set this up
    lines = []
    colors = []
    for edge in edges:
        lines.append([edge[0].loc, edge[1].loc])
        colors.append((1, 0, 0, 1))
    min_x = math.inf
    max_x = -math.inf
    min_y = math.inf
    max_y = -math.inf
    for point in points:
        if point.loc[0] < min_x:
            min_x = point.loc[0]
        if point.loc[0] > max_x:
            max_x = point.loc[0]
        if point.loc[1] < min_y:
            min_y = point.loc[1]
        if point.loc[1] > max_y:
            max_y = point.loc[1]
        for arc in point.arcs:
            lines.append([point.loc, arc.loc])
            colors.append((0, 1, 0, 1))
    for index in range(len(path) - 1):
        lines.append([path[index].loc, path[index + 1].loc])
        colors.append((0, 0, 1, 1))
    graph = mpl.LineCollection(lines, colors=colors, linewidths=2)
    fig, ax = plt.subplots()
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)
    ax.add_collection(graph)
    plt.show()

def main():
    """
    main function of the program
    """
    # points = [Vertex(1,1), Vertex(2,2), Vertex(3,0), Vertex(4,5), Vertex(5,3), Vertex(6,9)]
    # edges = [(points[0], points[1]), (points[0], points[2]), (points[2], points[1]), (points[3], points[4]), (points[3], points[5]), (points[4], points[5])]
    points, edges = generate_obstacles(210, "nlogn")
    start = Vertex(0, 0)
    end = random.choice(points)
    naive_visibility_graph(points, edges, start)
    points.append(start)
    shortest_path = dijkstras_alg(points, edges, start, end)
    test_arcs(points, edges, shortest_path)
    
main()

