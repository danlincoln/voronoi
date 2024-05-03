from itertools import product
from random import sample

import click
import numpy as np

from src.guibas_and_stolfi import delaunay
from src.dual import dfs_edges, calculate_dual

rng = np.random.default_rng()


@click.command()
@click.option("-h", "--height", default=1080, help="Image height (default 1080).")
@click.option("-w", "--width", default=1920, help="Image width (default 1920).")
@click.option(
    "-s",
    "--sites",
    "num_sites",
    default=10,
    help="Number of points to use for the diagram (default 10).",
)
@click.option(
    "--foreground",
    default="rgb(0, 0, 0)",
    help="Foreground cell color in RGB format (default black).",
)
@click.option(
    "--background",
    default="rgb(255, 255, 255)",
    help="Background cell color in RGB format (default white).",
)
@click.option(
    "--line",
    default="rgb(255, 255, 255)",
    help="Line color in RGB format (default white).",
)
@click.option("--line_weight", default=1, help="Line width (default 1).")
@click.argument("file", type=click.File("w"))
def diagram(
    height: int,
    width: int,
    num_sites: int,
    foreground: str,
    background: str,
    line: str,
    line_weight: int,
    file: str,
) -> None:
    """Calculate a Voronoi diagram for a number of random two-dimensional sites. Output
    the results as an SVG.

    FILE is the name of the file to write.
    """
    sites = sample(list(product(range(width), range(height))), k=num_sites)
    # Add four sites at extremities to make the image look nicer at the edges.
    sites.extend([(-5, -5), (-5, height + 5), (width + 5, height + 5), (width + 5, -5)])
    sites = sorted(sites)
    l, _ = delaunay(sites)
    calculate_dual(l)

    file.write(
        f'<svg height="{height}" width="{width}" '
        + 'xmlns="http://www.w3.org/2000/svg">\n'
    )
    file.write(
        f'\t<rect height="{height}" width="{width}" style="fill:{background}" />\n'
    )
    for e in dfs_edges(l):
        points = " ".join(f"{x},{y}" for x, y in e.data["vertices"])
        color = f"{foreground[:-1]}, {rng.random() / 2})"
        file.write(
            f'\t<polygon points="{points}" '
            + f'style="fill:{color}; '
            + f"stroke:{line}; "
            + f'stroke-width:{line_weight}" />\n'
        )
    file.write("</svg>")


if __name__ == "__main__":
    diagram()
