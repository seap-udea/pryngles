##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#

# PlanetaRY spanGLES                                             #
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
# Main contributors:                                             #
#   Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado         #
##################################################################
#!/bin/bash
. .pack/packrc
buildpath="src/$PACKNAME/.build"

##################################################################
# Preliminary
##################################################################
if [ "x$1" = "x" ];then 
    echo "You must provide at least one .ipynb file"
fi

GIT=0
if [ -d .git ];then
    GIT=1
fi
IFS_ORIGINAL=$IFS

##################################################################
# Routines
##################################################################
function newer()
{
    file1=$1;shift
    file2=$1;shift

    if [ ! -e $file1 ];then
	echo -1
	return
    fi

    if [ ! -e $file2 ];then
	echo 1
	return 
    fi

    dif=$(($(date -r $file1 +%s)-$(date -r $file2 +%s)))
    echo $dif
    return
}

##################################################################
# Conversion
##################################################################
forced=0
for notebook in $@
do
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Check notebook
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if [ ! -e $notebook ];then
	if [ $notebook == "forced" ];then
	    forced=1
	else
	    echo "Notebook $notebook does not exist. Skipping."
	fi
	continue
    fi

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Parse module name
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    devfile=$(basename $notebook)
    devdir=$(dirname $notebook)

    IFS="-"
    targetdir="src"
    for dir in $devfile
    do
	if [ -d $targetdir/$dir ];then
	    targetdir="$targetdir/$dir"
	else
	    filename=$dir
	fi
    done
    IFS=" "
    filename=$(echo $filename |awk -F'.' '{print $1}')

    if ! [[ $notebook == *"$PACKNAME-"* ]]
    then 
	target=$devdir/$filename.py
	test=$devdir/test-$filename.py
    else
	target=$targetdir/$filename.py
	test=$targetdir/tests/test-$filename.py
    fi

    echo "xconverting $notebook into $target ($test)..."

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Update source file
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if [ $(newer $notebook $target) -gt 0 ]
    then
	bash bin/convert.sh $notebook
    fi

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create source files
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #Module template file
    template=$filename
    if [ ! -e src/$template.temp ]
    then
	template=""
    fi
    (cat src/header.py;cat src/$template.temp) > $buildpath/$filename.temp

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Parse source file
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    source=$buildpath/$filename.src
    python bin/jupdev-parse.py $source
    
done

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Build modules
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFS=$IFS_ORIGINAL
for module in $(cat src/modules)
do
    echo "Building module $module..."
    srcpath=$buildpath/$module
    modpath=src/$PACKNAME/$module.py
    (cat $srcpath.temp;cat $srcpath.external;cat $srcpath.standalone) > $modpath

    #Add classes
    for class in $(cat $srcpath.classes)
    do
	echo -e "\tAdding class $class..."
	cat ${srcpath}_$class.class >> $modpath

	#Add methods of class
	for methods in $(ls $buildpath/*_$class.methods 2>/dev/null)
	do
	    echo -e "\t\tAdding methods from $methods..."
	    cat $methods >> $modpath
	done
    done
done

#Compile constants
constpath=src/$PACKNAME/consts.py
echo "Adding constants..."
for module in $(cat src/modules)
do
    for modconst in $(ls $buildpath/*_$module.consts 2>/dev/null)
    do
	echo -e "\tAdding constants from $modconst..."
	cat $modconst >> $constpath
    done
done
#for modulepath in dev/pryngles-*.ipynb
#filename=$(basename $modulepath)
#module=$(echo $filename |awk -F"." '{print $1}' |awk -F"pryngles-" '{print $2}')
