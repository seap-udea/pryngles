#!/usr/bin/env python
##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
from sys import argv
import os

if len(argv)<=1:
    raise AssertionError("You must provide a filename.")
elif not os.path.exists(argv[1]):
    raise ValueError(f"Filename provided '{argv[1]}' not found.")
filepath=argv[1]
filename=os.path.basename(filepath)
filetest=f"/tmp/test-{filename}"
filesrc=f"/tmp/{filename}"

#Create
os.system(f"cat header.py > {filetest}")
fp=open(filepath,"r")
fo=open(filesrc,"w")
ft=open(filetest,"a")

ft.write("""import unittest
from pryngles import *
class Test(unittest.TestCase):
""")

qtest=False
qcont=False
for line in fp:
    if "IN_JUPYTER_TEST" in line:
        qtest=True
    if qtest:
        if "def test_" in line:
            qcont=True
        if "class Test" in line:
            qcont=False
        if qcont:
            ft.write(line)
    else:
        fo.write(line)

    if "unittest.main" in line:
        qtest=False


ft.write("""
if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
""")
        
fp.close()
fo.close()
ft.close()
