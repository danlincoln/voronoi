"""Calculate a Delaunay triangulation using Guibas & Stolfi's divide-and-conquer
algorithm. This file implements the Delaunay triangulation algorithm and the authors'
edge algebra and quad-edge data structure.

Guibas, L., & Stolfi, J. (1985). Primitivies for the Manipulation of General
    Subdivisions and the Computation of Voronoi Diagrams. ACM Transactions on Graphics,
    Vol. 4 (No. 2), pp. 74-123. DOI: https://doi.org/10.1145/282918.282923
"""

import numpy as np


class Quad:
    """Represents a set of four related Edges: the Edge, its symmetric counterpart, its
    dual, and its dual's symmetric counterpart. Every Edge in a Quad is a rotation of
    the first Edge in the Quad.
    """

    def __init__(self) -> None:
        self.edges = []

    def __getitem__(self, n: int) -> "Edge":
        """Returns the quarter-Edge at `n`."""
        return self.edges[n]

    @staticmethod
    def make() -> "Edge":
        """Create and return a new subdivision (Guibas & Stolfi, 1985, p.96)."""
        quad = Quad()

        # A new subdivision is just a line (i.e., org != dest and next = self).
        edge = Edge(0, quad)
        edge.next = edge

        sym = Edge(2, quad)
        sym.next = sym

        # The dual edges form a loop (i.e., org = dest and next = sym).
        rot = Edge(1, quad)
        rotsym = Edge(3, quad)
        rot.next = rotsym
        rotsym.next = rot

        quad.edges = [edge, rot, sym, rotsym]

        return edge

    def __repr__(self) -> str:
        return f"Quad{self[0].__repr__()}"


class Edge:
    def __init__(self, rot: int, quad: Quad) -> None:
        self._rot = rot
        self._quad = quad
        self.next = None
        self.data = {"origin": None}

    # ################ Fundamental Methods ################
    @property
    def org(self) -> object:
        """Returns the Edge's origin data."""
        return self.data["origin"]

    @org.setter
    def org(self, d: object) -> None:
        """Sets the Edge's origin data."""
        self.data["origin"] = d

    def rot(self, n: int) -> "Edge":
        """Returns the edge obtained by rotating the given edge `n` times. Each Quad
        contains the four rotated versions of each edge. The given edge's rotation
        (`e._rot`) added to `n` modulo 4 is used to determine which of the four edges to
        request form the Quad.
        """
        n = (self._rot + n) % 4
        return self._quad[n]

    # ################## Derived Methods ##################
    @property
    def sym(self) -> "Edge":
        """Returns the given Edge's symmetrical pair. The origin of `e.sym` is the
        destination of `e` and vice-versa (Guibas & Stolfi, 1985, p.80).
        """
        return self.rot(2)

    @property
    def dest(self) -> object:
        """Returns the Edge's destination data."""
        return self.sym.org

    @dest.setter
    def dest(self, d: object) -> None:
        """Sets the Edge's destination data."""
        self.sym.org = d

    @property
    def left(self) -> object:
        """Returns the Edge's left face."""
        return self.rot(-1)

    @property
    def right(self) -> object:
        """Returns the Edge's right face."""
        return self.rot(1)
    
    @property
    def onext(self) -> "Edge":
        """Returns the next (counter-clockwise) Edge sharing the same left face and
        origin. Alias for `e.next` (Guibas & Stolfi, 1985, p.81).
        """
        return self.next

    @property
    def lnext(self) -> "Edge":
        """Returns the next (counter-clockwise) Edge sharing the same left face (Guibas
        & Stolfi, 1985, p.81).
        """
        return self.rot(-1).onext.rot(1)

    @property
    def rnext(self) -> "Edge":
        """Returns the next (counter-clockwise) Edge sharing the same right face.
        (Guibas & Stolfi, 1985, p.84).
        """
        return self.rot(1).onext.rot(-1)

    @property
    def dnext(self) -> "Edge":
        """Returns the next (counter-clockwise) Edge sharing the same right face and
        destination (Guibas & Stolfi, 1985, p.84).
        """
        return self.sym.onext.sym

    @property
    def rprev(self) -> "Edge":
        """Returns the previous (clockwise) Edge sharing the same right face (Guibas &
        Stolfi, 1985, p.84).
        """
        return self.sym.onext

    @property
    def oprev(self) -> "Edge":
        """Returns the previous (clockwise) Edge sharing the same right face and origin
        (Guibas & Stolfi, 1985, p.84).
        """
        return self.rot(1).onext.rot(1)

    @property
    def lprev(self) -> "Edge":
        """Returns the previous (clockwise) Edge sharing the same left face (Guibas &
        Stolfi, 1985, p.84).
        """
        return self.onext.sym

    @property
    def rprev(self) -> "Edge":
        """Returns the previous (clockwise) Edge sharing the same right face (Guibas &
        Stolfi, 1985, p.84).
        """
        return self.sym.onext

    @property
    def dprev(self) -> "Edge":
        """Returns the previous (clockwise) Edge sharing the same left face and
        destination (Guibas & Stolfi, 1985, p.84).
        """
        return self.rot(-1).onext.rot(-1)

    # ################ Supporting Methods #################
    def __repr__(self) -> str:
        return f"<{self.org}->{self.dest}>"

    def __hash__(self) -> int:
        """Use the hash of the parent quad in place of its child Edges' hashes."""
        return hash(self._quad)  # Treat the parent Quad as all the siblings' ID.

    def __eq__(self, other: "Edge") -> bool:
        """Edges and their symmetric counterparts are equivalent for identity purposes."""
        return (self.org == other.org and self.dest == other.dest) or (
            self.org == other.dest and self.dest == other.org
        )

    def delete(self) -> None:
        """Clean up references from this Edge and all its siblings, then delete the
        parent Quad.
        """
        quad = self._quad
        other_edges = self._quad.edges[: self._rot] + self._quad.edges[self._rot + 1 :]
        for edge in other_edges:
            edge.next = None
            edge._quad = None
        self.next = None
        self._quad = None
        del quad


def splice(a: Edge, b: Edge) -> None:
    """Splice two edge rings, achieving:

    1. If the two rings are distinct, combine them into one.
    2. If the two are the same ring, break it into separate pieces.

    (Guibas & Stolfi, 1985, p.96)
    """
    alpha = a.next.rot(1)
    beta = b.next.rot(1)
    alpha.next, beta.next = beta.next, alpha.next
    a.next, b.next = b.next, a.next


def connect(a: Edge, b: Edge) -> "Edge":
    """Add a new Edge `e` connecting the destination of `a` to the origin of `b` in such
    a way that `a.left` = `b.left` = `e.left` after the connection is complete (Guibas &
    Stolfi, 1985, p.103).
    """
    e = Quad.make()
    e.org = a.dest
    e.dest = b.org
    splice(e, a.lnext)
    splice(e.sym, b)
    return e


def delete_edge(e: Edge) -> None:
    """Disconnect Edge `e` from the rest of the structure. This may cause the rest of
    the structure to fall apart into separate components (Guibas & Stolfi, 1985, p.103).
    """
    splice(e, e.oprev)
    splice(e.sym, e.sym.oprev)
    e.delete()


def in_circle(a: tuple, b: tuple, c: tuple, d: tuple) -> bool:
    """Return True if and only if point `d` is interior to the region of the plane that
    is bounded by the oriented circle `abc` and lies to the left of it (Guibas & Stolfi,
    1985, p.106).
    """
    # If the points are co-circular, skip finding the determinant and return false. This
    # check counts distinct points to avoid false positives from floating point errors.
    if len(set([a, b, c, d])) == 4:
        m = np.array(
            [
                [a[0], a[1], pow(a[0], 2) + pow(a[1], 2), 1],
                [b[0], b[1], pow(b[0], 2) + pow(b[1], 2), 1],
                [c[0], c[1], pow(c[0], 2) + pow(c[1], 2), 1],
                [d[0], d[1], pow(d[0], 2) + pow(d[1], 2), 1],
            ]
        )
        return np.linalg.det(m) > 0
    else:
        return False


def ccw(a: tuple, b: tuple, c: tuple) -> bool:
    """Return true if the points `a`, `b`, and `c` form a counter-clockwise oriented
    triangle (Guibas & Stolfi, 1985, p.113).
    """
    m = np.array(
        [
            [a[0], a[1], 1],
            [b[0], b[1], 1],
            [c[0], c[1], 1],
        ]
    )
    return np.linalg.det(m) > 0


def right_of(x: tuple, e: Edge) -> bool:
    """Wraps `ccw()` to determine if `x` and `e` represent a clockwise triangle."""
    return ccw(x, e.dest, e.org)


def left_of(x: tuple, e: Edge) -> bool:
    """Wraps `ccw()` to determine if `x` and `e` represent a counter clockwise triangle."""
    return ccw(x, e.org, e.dest)


def delaunay(S: list[tuple], op_queue: list = None) -> tuple[Edge, Edge]:
    """Recursively compute the Delaunay triangulation for each half of a set of sites
    `S` before marrying the two half triangulations into the triangulation for the
    entire set.

    Args:
        S: The list of sites over which to calculate the Delaunay triangulation.
        op_queue: Optional list (queue) to contain edge add/remove operations over time.
            This queue can be used to visualize the algorithm's progress in finding the
            triangulation.

    Returns:
        A tuple of two edges representing the left and right outer directed edges of the
        Delaunay triangulation over the given sites.

    (Guibas & Stolfi, 1985, p.114)
    """
    has_queue = op_queue is not None
    if len(S) == 2:
        # S contains two sites in sorted order (s1, s2).
        # Create an edge from s1 to s2.
        a = Quad.make()
        a.org, a.dest = S
        if has_queue:
            op_queue.append(("add", (a.org, a.dest)))
        return a, a.sym
    if len(S) == 3:
        # S contains three sites in sorted order (s1, s2, s3).
        # Create edges a connecting s1 to s2 and b connecting s2 to s3.
        s1, s2, s3 = S
        a, b = Quad.make(), Quad.make()
        splice(a.sym, b)
        a.org = s1
        a.dest = b.org = s2
        b.dest = s3
        if has_queue:
            op_queue.append(("add", (a.org, a.dest)))
            op_queue.append(("add", (b.org, b.dest)))
        # Now close the triangle.
        if ccw(s1, s2, s3):
            c = connect(b, a)
            if has_queue:
                op_queue.append(("add", (c.org, c.dest)))
            return a, b.sym
        elif ccw(s1, s3, s2):
            c = connect(b, a)
            if has_queue:
                op_queue.append(("add", (c.org, c.dest)))
            return c.sym, c
        else:
            # The three points are co-linear.
            return a, b.sym
    else:
        # S contains four or more poinrts. Let L and R be the left and right halves of S
        half = len(S) // 2
        L, R = S[:half], S[half:]
        ldo, ldi = delaunay(L, op_queue)
        rdi, rdo = delaunay(R, op_queue)

        # Compute the lower common tangent of L and R.
        while True:
            if left_of(rdi.org, ldi):
                ldi = ldi.lnext
            elif right_of(ldi.org, rdi):
                rdi = rdi.rprev
            else:
                break

        # Create a first cross edge base_l from rdi.org to ldi.org.
        base_l = connect(rdi.sym, ldi)
        if has_queue:
            op_queue.append(("add", (base_l.org, base_l.dest)))
        if ldi.org == ldo.org:
            ldo = base_l.sym
        if rdi.org == rdo.org:
            rdo = base_l

        # This is the merge loop.
        while True:
            # Locate the first L point (l_cand.dest) to be encountered by the rising
            # bubble, and delete L edges out of base_l.dest that fail the circle test.
            l_cand = base_l.sym.onext
            if right_of(l_cand.dest, base_l):  # Is l_cand valid?
                while in_circle(
                    base_l.dest, base_l.org, l_cand.dest, l_cand.onext.dest
                ):
                    t = l_cand.onext
                    if has_queue:
                        op_queue.append(("remove", (l_cand.org, l_cand.dest)))
                    delete_edge(l_cand)
                    l_cand = t

            # Symmetrically, locate the first R points to be hit and delete R edges.
            r_cand = base_l.oprev
            if right_of(r_cand.dest, base_l):  # Is r_cand valid?
                while in_circle(
                    base_l.dest, base_l.org, r_cand.dest, r_cand.oprev.dest
                ):
                    t = r_cand.oprev
                    if has_queue:
                        op_queue.append(("remove", (r_cand.org, r_cand.dest)))
                    delete_edge(r_cand)
                    r_cand = t

            # If both l_cand and r_cand are invalid, then base_l is the upper common
            # tangent.
            if not right_of(l_cand.dest, base_l) and not right_of(r_cand.dest, base_l):
                break

            # The next cross edge is to be connected to either l_cand.dest or
            # r_cand.dest. If both are valid, choose the appropriate one using the
            # circle test.
            if not right_of(l_cand.dest, base_l) or (
                right_of(r_cand.dest, base_l)
                and in_circle(l_cand.dest, l_cand.org, r_cand.org, r_cand.dest)
            ):
                # Add cross edge base_l from r_cand.dest to base_l.dest.
                base_l = connect(r_cand, base_l.sym)
                if has_queue:
                    op_queue.append(("add", (base_l.org, base_l.dest)))
            else:
                base_l = connect(base_l.sym, l_cand.sym)
                if has_queue:
                    op_queue.append(("add", (base_l.org, base_l.dest)))

        return ldo, rdo
