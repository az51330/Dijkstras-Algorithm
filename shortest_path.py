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

def euclidean_distance(point1, point2):
    """
    Calculates the euclidean distance of two given points
    
    keyword arguments:
    point1 -- tuple of x and y
    point2 -- tuple of x and y

    returns the float distance of the two points
    """
    return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
