"""
    File:        shortest_path.py
    Author:      A.J. Zuckerman
    Course:      CS 307 - Computational Geometry
    Assignment:  Synthesis - Shortest Path in Polygons
    Description: Finds the shortest path from start to goal point
    in a polygon with obstacles using Dijkstra's algorithm, runs experiment of
    query time to show O(nlogn + k) time, builds structure brute force since not testing that
"""

import math

class Vertex:
    def __init__(self, x, y):
        self.loc = (x, y)
        # an arc is a connection of two visible vertices listed by self loc then connect loc
        self.arcs = []

    def __str__(self):
        return f"({self.loc[0]}, {self.loc[1]})"

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

def naive_visibility_graph(points, edges, start):
    """
    Creates a visibility graph by brute force in n^3 time

    keyword arguments:
    points -- a list of points in the visibility graph (points on obstacles)
    edges -- all edges of obstacles
    start -- the start vertex

    returns the start vertex which connects to the rest of the graph through its arcs
    """
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
    return start

def test_arcs(points, start):
    for point in points:
        print(f"Point: {point}")
        print("Arcs:")
        for arc in point.arcs:
            print(arc)
    print(f"Point: {start}")
    print("Arcs:")
    for arc in start.arcs:
        print(arc)

def main():
    """
    main function of the program
    """
    points = [Vertex(1,1), Vertex(2,2), Vertex(3,0), Vertex(4,5), Vertex(5,3), Vertex(6,9)]
    edges = [(points[0], points[1]), (points[0], points[2]), (points[2], points[1]), (points[3], points[4]), (points[3], points[5]), (points[4], points[5])]
    start = Vertex(0, 4)
    naive_visibility_graph(points, edges, start)
    test_arcs(points, start)
    
main()
