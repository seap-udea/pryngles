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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
import setuptools
from numpy.distutils.core import Extension, setup

"""

Note about setuptools:

NumPy distutils is not simply the old distutils package.  Actually it
internally tries to import setuptools and use it as much as possible.
We use numpy.distutils because we are distributing Pryngles along with
the fortran code pixx.

For a detailed explanation see eg:
https://het.as.utexas.edu/HET/Software/Numpy/f2py/distutils.html
https://stackoverflow.com/a/55358607

"""

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='pryngles',
    author="Jorge I. Zuluaga, Allard Veenstra, Jaime A. Alvarado, Mario Sucerquia",
    author_email="jorge.zuluaga@udea.edu.co",
    description="PlanetaRY spanGLES: general photometry of planets, rings and shadows",
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
    version='0.9.5',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # EXTENSIONS
    # ######################################################################
    ext_modules=[
        Extension(name="pryngles.pixx",
                  sources=["src/pryngles/pixx_sig.pyf",
                           "src/pryngles/pixx/pixx.f",
                           ]),
        Extension(name="pryngles.cpixx",
                  sources=["src/pryngles/cpixx/cpixx.c",
                           ]),
    ],
    
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
    install_requires=['rebound','scipy','ipython',
	              'matplotlib','tqdm','dill',
	              'spiceypy','cmasher','pandas','celluloid',
	              'sigfig','anytree','ipywidgets'],

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={"": ["*.c","data/*.*", "tests/*.*"]},
)
