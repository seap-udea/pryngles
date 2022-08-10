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

##################################################################
#VARIABLES
##################################################################
PACKDIR=.pack/
include $(PACKDIR)/packrc
DEVFILES=$(shell ls dev/$(PACKNAME)-*.ipynb)
BRANCH=$(shell bash .getbranch.sh)
VERSION=$(shell tail -n 1 .versions)
PUBLIC=../$(PACKNAME)-public/
MOD=None

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
	@-find . -name "*.pyc" -delete
	@-find . -name '__pycache__' -type d | xargs rm -fr

cleandist:
	@-rm -rf dist/
	@-rm -rf build/

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
	@bash convert.sh $(DEVFILES)

#Example: make release RELMODE=release VERSION=0.2.0.2 
release:
	@echo "Releasing a new version..."
	@bash release.sh $(RELMODE) $(VERSION)

install:
	@echo "Installing system dependencies..."
	@bash .pack/deps.sh .pack/deps_pack_$(SYSTEM).conf
	@echo "Installing locally..."
	@$(PIP) install -e .

test:convert
	@echo "Testing package..."
ifeq ($(MOD),None)
	@$(NOSETESTS)
else
	@$(NOSETESTS) src/pryngles/tests/test-$(MOD).py
endif

version:
	@pip show $(PACKNAME)

public:
	@echo "Updating public github repo..."
	@cp src/pryngles/*.py $(PUBLIC)/src/pryngles
	@cp src/pryngles/tests/*.py $(PUBLIC)/src/pryngles/tests
	@cp -rf src/pryngles/data $(PUBLIC)/src/pryngles
	@cp examples/pryngles-tutorial-quickstart.ipynb $(PUBLIC)/
	@cp examples/pryngles-tutorial-developers.ipynb $(PUBLIC)/
	@cp papers/bright-side/pryngles-paper-figures.ipynb examples/pryngles-examples-exploration.ipynb
	@cp examples/pryngles-examples-exploration.ipynb $(PUBLIC)/
	@cp README.md LICENSE $(PUBLIC)/
	@make -C $(PUBLIC) commit

prueba:
