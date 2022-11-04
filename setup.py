##################################################################
#                                                                #
# #####...#####...##..##..##..##...####...##......######...####..#
# ##..##..##..##...####...###.##..##......##......##......##.....#
# #####...#####.....##....##.###..##.###..##......####.....####..#
# ##......##..##....##....##..##..##..##..##......##..........##.#
# ##......##..##....##....##..##...####...######..######...####..#
# ...............................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
import setuptools
from numpy.distutils.core import setup, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

pixx = Extension(name="pryngles.pixx",
                        sources=["src/pryngles/pixx_sig.pyf",
                                 "src/pryngles/pixx/reflection.f",
                                 "src/pryngles/pixx/rdfous_planet.f",
                                 "src/pryngles/pixx/rdfous_ring.f",
                                 "src/pryngles/pixx/bracks.f",
                                 "src/pryngles/pixx/spline.f",
                                 "src/pryngles/pixx/splint.f",
                                 ])

setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='pryngles',
    author="Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado",
    author_email="jorge.zuluaga@udea.edu.co",
    description="PlanetaRY spanGLES: the bright-side of the light-curve of (ringed) exoplanets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/pryngles",
    keywords='astronomy exoplanets planetary-rings',
    license='MIT',

    # ######################################################################
    # CLASSIFIER
    # ######################################################################
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    version='0.9.0',

    # ######################################################################
    # FILES
    # ######################################################################
#    packages=["src"],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # EXTENSIONS
    # ######################################################################
    ext_modules=[pixx],
    
    # ######################################################################
    # ENTRY POINTS
    # ######################################################################
    entry_points={
        'console_scripts': ['install=pryngles.install:main'],
    },

    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='nose.collector',
    tests_require=['nose'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=[
        'rebound', 'scipy', 'ipython', 'matplotlib', 'tqdm',
        'dill', 'spiceypy', 'cmasher','pandas','celluloid',
        'sigfig','anytree'
    ],

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={"": ["data/*.*", "tests/*.*"]},
)
