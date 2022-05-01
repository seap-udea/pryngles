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

if [ "x$1" = "x" ];then 
    echo "You must provide at least one .ipynb file"
fi
GIT=0
if [ -d .git ];then
    GIT=1
fi

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

function convert()
{
    notebook=$1;shift
    target=$1;shift
    test=$1;shift
    
    echo -e "\tConverting from ipynb $notebook to python $target..."
    jupyter nbconvert --to python $notebook --stdout 2> /dev/null | grep -v "# In" | cat -s > /tmp/convert.py 

    echo -e "\tProcessing magic commands..."
    sed -ie "s/get_ipython().magic('timeit\(.*\))$/get_ipython().magic('timeit\1,scope=globals())/" /tmp/convert.py

    echo -e "\tTriming end..."
    nlines=$(cat -n /tmp/convert.py | grep -e "--End--" | cut -f 1 )
    if [ "x$nlines" = "x" ];then
	nlines=$(cat /tmp/convert.py|wc -l)
    else
	((nlines--))
    fi

    echo -e "\tTriming test..."
    ntlines=$(cat -n /tmp/convert.py | grep -e "--Test--" | cut -f 1 )
    if [ "x$nlines" = "x" ];then
	ntlines=$(cat /tmp/convert.py|wc -l)
    else
	((ntlines--))
    fi

    if [ $ntlines -eq -1 ];then
	tlines=$nlines
    else
	tlines=$ntlines
    fi
    
    echo -e "\tAdding header..."
    (cat header.py;head -n $tlines /tmp/convert.py) > $target

    if [ $ntlines -gt 0 ];then
	alines=$(cat /tmp/convert.py|wc -l)
	elines=$((alines-ntlines))
	rlines=$((nlines-ntlines))
	(cat header.py;tail -n $elines /tmp/convert.py | head -n $rlines) > $test
	echo -e "\tCreating test file $test..."
    fi
	
}

for notebook in $@
do
    if [ ! -e $notebook ];then 
	echo "Notebook $notebook does not exist. Skipping."
	continue
    fi

    devfile=$(basename $notebook)
    devdir=$(dirname $notebook)

    # Parse script name
    IFS="-"
    targetdir="src/"
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

    # Check if notebook is more recent than target file
    if [ $1 != "force" ];then 
	if [ $(newer $notebook $target) -lt 0 ];then continue;fi
    fi

    echo "Analysing file $devfile:"
    if [ $GIT -gt 0 ];then
	git add -f $notebook
    fi

    echo -e "\tDirectory: $targetdir"
    echo -e "\tFilename: $filename"
    echo -e "\tTarget object: $target"

    convert $notebook $target $test

    if [[ $notebook == *"$PACKNAME-"* ]]
    then
	if [ $GIT -gt 0 ];then
	    git add -f $target
	fi
    fi
done
echo "Completed."
