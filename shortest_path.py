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
                    if orient(edge[0].loc, edge[1].loc, point.loc) == -1 and orient(edge[0].loc, edge[1].loc, other.loc) == 1:
                        if orient(edge[0].loc, edge[1].loc, point.loc) == 1 and orient(edge[0].loc, edge[1].loc, other.loc) == -1:
                            arc = False
                if arc:
                    point.arc.append(other)
    return start

def main():
    """
    main function of the program
    """
    print(Vertex(1,1))

main()
