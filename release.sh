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
#!/bin/bash
. .pack/packrc
cd pryngles

version=$1
setversion=$(grep "version=" setup.py |awk -F"'" '{print $2}')

##################################################################
# Latest version
##################################################################
if [ "x$version" = "x" ]
then
    version=$(tail -n 1 .versions)
    echo "Latest version: $version"
    echo "Version in setup file: $setversion"
    exit 1
fi

if [ "$version" = "$setversion" ]
then
    echo "Version provided ($version) coincide with version in setup.py file ($setversion). It must be different."
    exit 1
fi

echo "Releasing version $version of the package..."

##################################################################
# Update setup.py file
##################################################################
sed -i.bak "s/version=\'[0-9\.]*\'/version='$version'/gi" setup.py 

##################################################################
# Update files
##################################################################
echo "Updating package files..."
cp -rf ../README.md .
rm -rf src/pryngles/*
cp -rf *.py src/pryngles/
cp -rf data src/pryngles/
rm -rf src/pryngles/setup.py

##################################################################
# Remove previous versions
##################################################################
echo "Removing previous version..."
rm -rf dist/*

##################################################################
# Build package
##################################################################
echo "Building packages..."
python -m build

##################################################################
# Uploading the package
##################################################################
echo "Uploading to PyPI (use __token__ as username and pypi-<token> as password..."
python -m twine upload --repository testpypi dist/* --verbose

##################################################################
# Report version
##################################################################
echo $version >> .versions
cd -
