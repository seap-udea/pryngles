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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
"""This file parses the python code of a module searching for "inline"
test code.

Inline test code is code inserted into the source Jupyter notebooks of
module files using the diective #@test:plot...#@end:test.

Arguments:

  Name of the file to convert.

"""
from sys import argv
import os

if len(argv)<=1:
    raise AssertionError("You must provide a filename.")
elif not os.path.exists(argv[1]):
    raise ValueError(f"Filename provided '{argv[1]}' not found.")

for filepath in argv[1:]:
    filename=os.path.basename(filepath)
    if "ipynb" not in filename:
        raise ValueError(f"You must provide a Python Notebook with extension .ipynb ('{filename}' provided).")
    filedir=os.path.dirname(filepath)
    # Parts of the name
    filebase=filename.split(".")[0]
    filesrc=f"/tmp/convert-{filebase}.py"
    filetar=f"/tmp/{filebase}.py"

    #Target file
    filetest=f"{filedir}/test-{filebase}.py"

    if os.path.exists(filetest):
        os.system(f"rm {filetest}")
    if os.path.exists(filesrc):
        os.system(f"rm {filesrc}")

    print(f"Extracting test code from {filepath} and storing it into {filetest}...")

    #Converting notebook
    print(f"\tConverting notebook into python file...")
    os.system(f"jupyter nbconvert --to python {filepath} --stdout 2> /dev/null | grep -v '# In' | cat -s > {filesrc}")

    #Create
    print(f"\tGetting the header...")
    os.system(f"cat src/pryngles/header.py > {filetest}")

    fp=open(filesrc,"r")
    fo=open(filetar,"w")
    ft=open(filetest,"a")

    print(f"\tWritting the file (the remaining code will be stored at {filetar})...")
    ft.write("""import unittest
from pryngles import *
class Test(unittest.TestCase):
""")

    qtest=False
    qfulltest=False
    qcont=False
    qfullcont=False
    ntest=0
    for line in fp:
        if "#@test:" in line:
            ntest+=1
            qtest=True

        if "#@fulltest" in line:
            qfulltest=True
            ntest=1

        if qtest:
            if ("def test_" in line):
                qcont=True
            if ("class Test" in line):
                qcont=False
            if qcont:
                ft.write("\t"+line)
        elif qfulltest:
            if ("class Test" in line) or ("@fulltest" in line) or ("@end:fulltest" in line):
                qfullcont=False
            else:
                qfullcont=True
            if qfullcont:
                ft.write(line)
        else:
            fo.write(line)

        if "@end:test" in line:
            qtest=False

    if not qfulltest:
        ft.write("""\
if __name__=="__main__":
   unittest.main(argv=['first-arg-is-ignored'],exit=False)
""")

    fp.close()
    fo.close()
    ft.close()

    print(f"Found {ntest} tests in source file...")

