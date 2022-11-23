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
function extract()
{
    file=$1;shift
    initial_string=$1;shift
    final_string=$1;shift
    
    #Number of lines
    nlines=$(wc -l $file |awk '{print $1}')
    
    #Get dependencies
    ndep_ini=$(grep -Hn $initial_string $file |head -n 1 |awk -F':' '{print $2}')
    
    #Get rest of file
    ndep_tail=$((nlines-ndep_ini+1))
    ndep_up=$(tail -n $ndep_tail $file |grep -Hn $final_string |head -n 1 |awk -F':' '{print $2}')
    
    #Get the lines of install
    dependencies_string=$(tail -n $ndep_tail $file |head -n $ndep_up)

    dependencies_list=()
    IFS=","
    j=0
    for dep in $(echo $dependencies_string | cut -d "[" -f2 | cut -d "]" -f1 |sed s/\'//gi)
    do
	package=$(echo $dep |sed -e 's/^[[:space:]]*//' |sed -e 's/^[[:space:]]*//' |sed 's/ *$//g')
	dependencies_list[$j]=$package
	((j++))
    done
    echo ${dependencies_list[@]}
}

####################################################################
# SEARCH FOR PYTHON
####################################################################
if python3 --version &> /dev/null
then
    PYTHON=python3
    PIP=pip3
elif python --version &> /dev/null
then
    PYTHON=python
    PIP=pip
else
    echo "No python installed"
    exit 1
fi

if nosetests3 --version &> /dev/null
then
    NOSETESTS=nosetests3
elif nosetests --version &> /dev/null
then
    NOSETESTS=nosetests
else
    $PYTHON -m pip install nose
fi

####################################################################
# INSTALL DEPENDENCIES
####################################################################
deps_setup=$(extract setup.py "install_requires" "]")
deps_pyproj=$(extract pyproject.toml "requires" "]")

$PYTHON -m pip install -q $deps_setup $deps_pyproj

####################################################################
# CREATE PACKAGE CONFIGURATION
####################################################################
echo "
##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES:                                            #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
#Pthon binaries. In Linux is preferrable python3, nosetests3
PYTHON=$PYTHON
NOSETESTS=$NOSETESTS
PIP=$PIP

#Directories
STOREDIR=.store

#Log files
LOGFILE=.pack/log

#System
SYSTEM=macosx
#SYSTEM=linux

#Other
DEVDIR=dev/
PACKNAME=pryngles

#Default mode of release
RELMODE=test
" > .pack/packrc
