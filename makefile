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

##################################################################
#VARIABLES
##################################################################
SHELL:=/bin/bash

#Package information
PACKDIR=.pack/
include $(PACKDIR)/packrc

#Github specifics
BRANCH=$(shell bash .getbranch.sh)
VERSION=$(shell tail -n 1 .versions)

#JupDev files
DEVFILES=$(shell ls dev/$(PACKNAME)-*.ipynb)

#Public repo
PUBLIC=../$(PACKNAME)-public/

#Module to convert
MOD=None

#Modules
MODULES=$(shell cat src/modules)

#Enforce conversion.  Use make convert ENFORCE=forced
ENFORCE=

show:
	@echo "Development files:" $(DEVFILES)
	@echo $(BRANCH)

##################################################################
#BRANCHING
##################################################################
branch:
	git checkout -b dev

make rmbranch:
	git branch -d dev

devbranch:
	git checkout dev

master:
	git checkout master

##################################################################
#BASIC RULES
##################################################################
deps:
	@echo "Checking dependencies for $(SYSTEM)..."
	@bash $(PACKDIR)/deps.sh $(PACKDIR)/deps_pack_$(SYSTEM).conf 

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

cleandist:
	@-rm -rf dist/
	@-rm -rf build/

cleanpack:
	@-cp src/$(PACKNAME)/version.py tmp/
	@-rm -rf src/$(PACKNAME)/*.py
	@-cp tmp/version.py src/$(PACKNAME)/version.py
	@-rm -rf src/$(PACKNAME)/tests/*.py
	@-rm -rf src/$(PACKNAME)/.build/*

##################################################################
#GIT
##################################################################
addall:cleanall
	@echo "Adding..."
	@-git add -A .

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin $(BRANCH)

pull:
	@echo "Pulling new files..."
	@-git reset --hard HEAD
	@-git pull origin $(BRANCH)

##################################################################
#PACKAGE RULES
##################################################################
pack:
	@echo "Packing data..."
	@bash $(PACKDIR)/pack.sh

unpack:
	@echo "Unpacking data..."
	@bash $(PACKDIR)/pack.sh unpack

convert:
	@echo "Converting iPython Notebooks $(DEVFILES)..."
	@bash bin/convert.sh $(ENFORCE) $(DEVFILES)

xconvert:convert
	@echo "Converting iPython Notebooks $(DEVFILES)..."
	@bash bin/xconvert.sh $(DEVFILES)

#Example: make release RELMODE=release VERSION=0.2.0.2 
release:
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)

install:
	@echo "Installing system dependencies..."
	@bash .pack/deps.sh .pack/deps_pack_$(SYSTEM).conf
	@echo "Installing locally..."
	@$(PIP) install -e .

pipinstall:
	@$(PIP) install -e .

import:
	@$(PYTHON) -c "from pryngles import *;print(version)"

test:import
ifeq ($(MOD),None)
	@echo "Testing all modules..."
	@$(NOSETESTS) 2> >(tee -a /tmp/$(PACKNAME)-test-errors.log >&2)
else
	@echo "Testing module(s) $(MOD)..."
	@for mod in $(shell echo $(MOD) | sed 's/,/ /'); do echo "Testing $$mod";$(NOSETESTS) src/pryngles/tests/test-$${mod}.py 2> >(tee /tmp/$(PACKNAME)-test-errors-$${mod}.log >&2);done
endif

testdet:
	@echo "Testing module by module..."
	@for module in $(MODULES);do make test MOD=$$module;done

version:
	@pip show $(PACKNAME)

public:
	@echo "Updating public github repo..."
	@cp src/pryngles/*.py $(PUBLIC)/src/pryngles
	@cp src/pryngles/tests/*.py $(PUBLIC)/src/pryngles/tests
	@cp -rf src/pryngles/data $(PUBLIC)/src/pryngles
	@cp examples/pryngles-tutorial-quickstart.ipynb $(PUBLIC)/
	@cp examples/pryngles-dev*-tutorial.ipynb $(PUBLIC)/
	@cp papers/bright-side/pryngles-paper-figures.ipynb examples/pryngles-examples-exploration.ipynb
	@cp examples/pryngles-examples-exploration.ipynb $(PUBLIC)/
	@cp README.md LICENSE WHATSNEW.md $(PUBLIC)/
	@make -C $(PUBLIC) commit
