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

function convert()
{
    notebook=$1;shift
    target=$1;shift
    test=$1;shift
    filename=$1;shift
    
    echo -e "\tConverting from ipynb $notebook to python $target..."
    jupyter nbconvert --to python $notebook --stdout 2> /dev/null | grep -v "# In" | cat -s > /tmp/convert.py 

    echo -e "\tProcessing magic commands..."
    sed -ie "s/get_ipython().magic('timeit\(.*\))$/get_ipython().magic('timeit\1,scope=globals())/" /tmp/convert.py

    echo -e "\tTriming end..."
    nlines=$(cat -n /tmp/convert.py | grep -e "#@end:module" | cut -f 1 )
    if [ "x$nlines" = "x" ];then
	nlines=$(cat /tmp/convert.py|wc -l)
    else
	((nlines--))
    fi

    #Deprecated
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
    
    #echo -e "\tAdding header..."
    if [ ! -e src/$filename.temp ]
    then
	filename=""
    fi
    echo -e "\tUsing as template src/$filename.temp"
    
    (cat src/header.py;cat src/$filename.temp;head -n $tlines /tmp/convert.py) > $target
    #(head -n $tlines /tmp/convert.py) > $target

    if [ $ntlines -gt 0 ];then
	alines=$(cat /tmp/convert.py|wc -l)
	elines=$((alines-ntlines))
	rlines=$((nlines-ntlines))
	(cat src/header.py;tail -n $elines /tmp/convert.py | head -n $rlines) > $test
	#(tail -n $elines /tmp/convert.py | head -n $rlines) > $test
	echo -e "\tCreating test file $test..."
    fi
    
}

##################################################################
# Convert
##################################################################
forced=0
for notebook in $@
do
    if [ ! -e $notebook ];then
	if [ $notebook == "forced" ];then
	    forced=1
	else
	    echo "Notebook $notebook does not exist. Skipping."
	fi
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
    if [ $forced -lt 1 ];then 
	if [ $(newer $notebook $target) -lt 0 ];then continue;fi
    fi

    echo "Analysing file $devfile:"
    if [ $GIT -gt 0 ];then
	git add -f $notebook
    fi

    echo -e "\tDirectory: $targetdir"
    echo -e "\tFilename: $filename"
    echo -e "\tTarget object: $target"
    filebase=$filename
    
    convert $notebook $target $test $filename

    # Parsing inline test code
    echo -e "\tParsing python file $target"
    ntests=$($PYTHON bin/test-parse.py $target)
    if [ $ntests -gt 0 ];then
	echo -e "\tParsing tests from inline code"
	cp -rf /tmp/test-$filebase.py $targetdir/tests/
	cp -rf /tmp/$filebase.py $target
    else
	echo -e "\tNo inline test code"
    fi

    #Save source file
    source=src/$PACKNAME/.build/$filebase.src
    cp -rf $target $source    

    # Add new package file into git repository
    if [[ $notebook == *"$PACKNAME-"* ]]
    then
	if [ $GIT -gt 0 ];then
	    git add -f $target
	    git add -f $targetdir/tests/test-$filebase.py
	fi
    fi

done
echo "Completed."
