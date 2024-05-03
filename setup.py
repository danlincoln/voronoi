from setuptools import setup, find_packages

setup(
    name="voronoi",
    version="1.0",
    py_modules=["voronoi"],
    packages=find_packages(include=["src"]),
    install_requires=[
        "Click",
        "Numpy"
    ],
    entry_points={
        'console_scripts': [
            'voronoi = voronoi:diagram',
        ],
    },
    license="MIT",
)
