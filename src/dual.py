"""Given a Delaunay triangulation, calculate its dual Voronoi diagram.

This file contains functions to calculate a Voronoi diagram's vertices and faces from
the output of Guibas and Stolfi's Delaunay triangulation algorithm.
"""

from src.guibas_and_stolfi import Edge
from math import acos, inf, sqrt, sin

import numpy as np


def calculate_dual(origin: Edge) -> None:
    """Given an `origin` edge of a Delaunay triangulation, calculate its dual (the
    Voronoi diagram) and add vertex data for both the Delaunay and Voronoi faces. The
    edges are updated in place.
    """

    # Calculate the outer edge ring.
    edge = origin
    outer_edges = [edge]
    while (edge := edge.rnext) != origin:
        outer_edges.append(edge)
    # The dual edges for the outer ring edges "start" at infinity—there are no points
    # beyond this outer ring to triangulate.
    data = {"origin": (inf, inf), "vertices": [e.org for e in outer_edges]}
    for e in outer_edges:
        e.right.data = data

    # Calculate all the interior triangles. For each face, check if its left and right
    # dual edges have been initialized. If not, add the vertices from the three
    # Delaunay edges that surround it, then calculate the origin using the circumcenter
    # of the triangle formed by those vertices.
    for edge in dfs_edges(origin):
        if edge.left.org is None:
            edge.left.data["vertices"] = [
                edge.org,
                edge.lnext.org,
                edge.lnext.lnext.org,
            ]
            edge.left.org = circumcenter(*edge.left.data["vertices"])

        if edge.right.org is None:
            edge.right.data["vertices"] = [
                edge.org,
                edge.rnext.org,
                edge.rnext.rnext.org,
            ]
            edge.right.org = circumcenter(*edge.right.data["vertices"])

    # Calculate a reasonable destination for dual edges that go to infinity. This makes
    # the resulting Voronoi diagram go all the way to the edges of the canvas. Otherwise
    # the outer edges might have gaps.
    for e in outer_edges:
        orthagonal_vector = (np.array(e.org) - np.array(e.dest))[::-1]
        orthagonal_vector *= np.array([-1, 1])
        dest = np.array(e.left.org) + orthagonal_vector
        data = {"origin": tuple(dest), "vertices": e.right.data["vertices"]}
        e.right.data = data

    # Add dual vertices to every primal edge, which allows easily getting the Voronoi
    # polygons for the visualization.
    for edge in dfs_edges(origin):
        vertices = [edge.left.org]
        while (edge := edge.onext).left.org not in vertices:
            vertices.append(edge.left.org)
        edge.data["vertices"] = vertices


def dfs_edges(e: Edge) -> set[Edge]:
    """Given a starting edge `e`, perform a depth-first search to find all connected
    edges in the graph.
    """
    stack: list[Edge] = list()
    discovered = set()
    stack.append(e)
    try:
        while e := stack.pop():
            if e not in discovered:
                discovered.add(e)
                stack.extend([e.rnext, e.onext, e.dnext, e.lnext])
    except IndexError:
        return discovered


def circumcenter(a: tuple, b: tuple, c: tuple) -> tuple:
    """Calculate and return the circumcenter point for a triangle represented by the
    vertices `a`, `b`, and `c`.
    """
    angles = (
        calc_angle(a, b, c),
        calc_angle(b, c, a),
        calc_angle(c, a, b),
    )
    xs, ys = [tuple(zip(z, angles)) for z in zip(a, b, c)]

    return (
        sum([point * sin(2 * angle) for point, angle in xs])
        / sum([sin(2 * angle) for angle in angles]),
        sum([point * sin(2 * angle) for point, angle in ys])
        / sum([sin(2 * angle) for angle in angles]),
    )


def calc_angle(a: tuple, b: tuple, c: tuple) -> float:
    """Given three points `a`, `b`, and `c` which form a triangle, calculate and return
    the angle at point `a` using the law of cosines.
    """
    A, B, C = dist(b, c), dist(a, c), dist(a, b)
    return acos((pow(B, 2) + pow(C, 2) - pow(A, 2)) / (2 * B * C))


def dist(a: tuple, b: tuple) -> float:
    """Calculate and return the distance between two points `a` and `b`."""
    return sqrt(pow(b[0] - a[0], 2) + pow(b[1] - a[1], 2))
