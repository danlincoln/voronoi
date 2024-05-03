# Voronoi

This program calculates a Delaunay triangulation and its Voronoi dual for arbitrary
points in two dimensions using the triangulation algorithm and quad-edge data structure
described by Leonidas Guibas and Jorge Stolfi in their 1985 paper
[*Primitivies for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams*](https://doi.org/10.1145/282918.282923).
It outputs the resulting diagram in SVG format.

## Usage

### Installation

1. Set up a Python virtual environment in the project directory.

        python -m venv venv

2. Activate the virtual environment.

        source venv/bin/activate

3. Install the module.

        pip install .

### Running

After the installation steps above, the application will run in your virtual environment
by calling `voronoi`.

There are several command-line options which can customize the output file.

        voronoi --help

        Usage: voronoi [OPTIONS] FILE

        Calculate a Voronoi diagram for a number of random two-dimensional sites.
        Output the results as an SVG.

        FILE is the name of the file to write.

        Options:
        -h, --height INTEGER   Image height (default 1080).
        -w, --width INTEGER    Image width (default 1920).
        -s, --sites INTEGER    Number of points to use for the diagram (default 10).
        --foreground TEXT      Foreground cell color in RGB format (default black).
        --background TEXT      Background cell color in RGB format (default white).
        --line TEXT            Line color in RGB format (default white).
        --line_weight INTEGER  Line width (default 1).
        --help                 Show this message and exit.

#### Example

For a green/purple diagram with dark purple lines, run:

        voronoi -h 200 -w 200 -s 50 --foreground "rgb(71, 28, 75)" --background "rgb(41, 98, 101)" --line "rgb(47, 34, 66)" output.svg

![A Voronoi diagram with purple cells on top of a green background.](img/output.svg)

## Formulas

The fomulas used in [src/dual.py](voronoi/dual.py) appear in mathematical notation below
because they helped my understanding. The formulas used in
[src/guibas_and_stolfi.py](voronoi/guibas_and_stolfi.py) appear in their paper: `in_circle`
on page 106 and `ccw` on page 113.

### Angle

`calc_angle` uses the law of cosines to find the interior angle of a point given two
others that form a triangle. Let $a$, $b$, and $c$ be three points $\in \mathbb{R}^2$
that form a triangle. Let $\alpha$ be the interior angle at point $a$.

$$
\alpha=\arccos\left(\frac{b^2+c^2-a^2}{2bc}\right)
$$

### Circumcenter

`circumcenter` calculates the circumcenter of a triangle formed by three points $a$,
$b$, and $c$ with interior angles $\alpha$, $\beta$, and $\gamma$ respectively
(calculated by `angle` above). Let $x$ be the triple formed by the three coordinates'
x-values, and $y$ be the triple formed by their y-coordinates. Then, the circumcenter
point for the triangle can be calculated by $O(x,y)$:

$$
O(x,y)=\left(
\frac{x_a\sin{2\alpha}+x_b\sin{2\beta}+x_c\sin{2\gamma}}
{\sin{2\alpha}+\sin{2\beta}+\sin{2\gamma}}
,
\frac{y_a\sin{2\alpha}+y_b\sin{2\beta}+y_c\sin{2\gamma}}
{\sin{2\alpha}+\sin{2\beta}+\sin{2\gamma}}
\right)
$$
