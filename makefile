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
PYTHON=python
PIP=pip
NOSE=nosetests
COMMIT="Update"
BRANCH=$(shell git branch |grep "*" |cut -f 2 -d " ")

##################################################################
#BASIC RULES
##################################################################
clean:cleancrap

cleanall:cleancrap cleanout cleandist

#=========================
#Clean
#=========================
cleancrap:
	@echo "Cleaning crap..."
	@-find . -name "*~" -delete
	@-find . -name "#*#" -delete
	@-find . -name "#*" -delete
	@-find . -name ".#*" -delete
	@-find . -name ".#*#" -delete
	@-find . -name ".DS_Store" -delete
	@-find . -name "Icon*" -delete
	@-find . -name "*.egg-info*" -type d | xargs rm -fr

cleanout:
	@echo "Cleaning all compiled objects..."
	@-find . -name "*.o" -delete
	@-find . -name "*.opp" -delete
	@-find . -name "*.gcno" -delete
	@-find . -name "*.gcda" -delete
	@-find . -name "*.gcov" -delete
	@-find . -name "*.info" -delete
	@-find . -name "*.out" -delete
	@-find . -name "*.tout" -delete
	@-find . -name "*.so" -delete
	@-find . -name '__pycache__' -type d | xargs rm -fr
	@-find . -name '.ipynb_checkpoints' -type d | xargs rm -fr

cleandist:
	@-rm -rf dist/
	@-rm -rf build/

##################################################################
#TEST
##################################################################
gentest:
	@echo "Generating tests..."
	@-$(PYTHON) bin/test-parse.py tests/*.ipynb

testall:
	@-$(NOSE) --verbosity=2 -x src/pryngles/tests

##################################################################
#GIT
##################################################################
push:
	@echo "Commiting..."
	git commit -am $(COMMIT)
	git push origin $(BRANCH)

pull:
	@echo "Pulling new files..."
	@-git pull

reset:
	@echo "Reset and pulling..."
	@-git --reset hard
	@-git pull

##################################################################
#PACKAGE RULES
##################################################################
#Example: make release RELMODE=release VERSION=0.2.0.2 
release:
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)

install:
	@$(PIP) install -e .

import:
	@$(PYTHON) -c "from pryngles import *;print(version)"
