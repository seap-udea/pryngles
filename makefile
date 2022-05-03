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

show:
	@echo "Development files:" $(DEVFILES)
	@echo $(BRANCH)

##################################################################
#BASIC RULES
##################################################################
deps:
	@echo "Checking dependencies for $(SYSTEM)..."
	@bash $(PACKDIR)/deps.sh $(PACKDIR)/deps_pack_$(SYSTEM).conf 

clean:cleancrap

cleanall:cleancrap cleanout

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

test:
	@echo "Testing package..."
	@$(NOSETESTS) $(PACKNAME)

version:
	@pip show $(PACKNAME)
