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
    dependencies=$(tail -n $ndep_tail $file |head -n $ndep_up)
    
    echo $dependencies
}

deps_setup=$(extract setup.py "install_requires" "]")
deps_pyproj=$(extract pyproject.toml "requires" "]")

echo $deps_setup > /tmp/deps.py
echo "print(install_requires)" >> /tmp/deps.py

$deps_name=$(echo $deps_setup | cut -d "[" -f2 | cut -d "]" -f1 |sed s/^[[:space:]]*// |xargs)

echo $deps_name

# for i in $(echo $deps_setup | cut -d "[" -f2 | cut -d "]" -f1)
# do
#     echo $i
#     #echo $i |sed s/\'//gi |xargs
# done
