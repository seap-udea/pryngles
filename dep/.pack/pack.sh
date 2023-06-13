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
#!/bin/bash
. .pack/packrc
basedir=$(pwd)

pack=$1
if [ "x$pack" == "x" ];then pack="pack";fi
confile="pack.conf"

GIT=0
if [ -d .git ];then
    GIT=1
fi

if [ $pack = "pack" ];then
    echo "Packing..."
    find $STOREDIR -name "*--*" -type d | xargs rm -rf 
    for file in $(cat $STOREDIR/$confile |grep -v "#")
    do
	echo -e "\tfile $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	uname=$(echo $dname |sed -e s/\\//_/)
	sdir="$STOREDIR/$uname--$fname"
	mkdir -p "$sdir"
	cd $sdir
	split -b 2000k $basedir/$file $fname-
	cd - &> /dev/null
	if [ $GIT -gt 0 ];then git add -f "$STOREDIR/$uname--$fname/";fi
    done
    if [ $GIT -gt 0 ];then 
	find $STOREDIR -name "*--*" -type d | xargs git add -f
    fi
else
    echo "Unpacking..."
    for file in $(cat $STOREDIR/$confile |grep -v "#")
    do
	echo -e "\tUnpacking $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	uname=$(echo $dname |sed -e s/\\//_/)
	sdir="$STOREDIR/$uname--$fname"
	cat "$sdir"/$fname-* > $dname/$fname
    done
fi
