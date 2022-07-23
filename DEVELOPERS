Guide to developers
===================

In this guide you will find common procedures for developers.

Version management
------------------

- Versions are in .versions file.

- Version 0: first release and design.  Example: 0.2.1

- Version 1: this will be the release when modularity of the package
  be improved.

- Version advance according to major changes.  Version 0.2.x are a
  family of the main version 0.2.

- All test versions have number.  Thus for instance: if we are in
  version 0.2.1, and want to release a test version it will be 0.2.1.1

- New version number scheme: 0.x.y (x-major, y-minor release), 0.x.y.z
  (z test version).


Installing locally pryngles
---------------------------

To install pryngles from the local directory:

```
python -m pip install -e .
```


Release a new version
---------------------

To release a new version of the package to the PyPI repo use:

```
make release RELMODE=release VERSION=x.y.z
```

where `x.y.z` is the new version.  The latest version can be obtained with:

```
make version
```

Release a test version
----------------------

To release a test version use:

```
make release RELMODE=test VERSION=x.y.z.w
```

Create a virtual environment to test Pryngles
---------------------------------------------

Virtual environments are the ideal place to test if a given version of
Pryngles is ready for being released.

You will need the package:

```
python -m pip install virtualenvwrapper
```

Locate the `virtualenvwrapper.sh` script:

```
which virtualenvwrapper.sh
```

For instance the location is `/usr/local/bin/virtualenvwrapper.sh`

Add the following line to .bashrc or .profile:

```
source /usr/local/bin/virtualenvwrapper.sh
```

To create a virtual environment use:

```
mkvirtualenv pryngles-lab
```

Now you need to create a kernelspec for this virtual environment:

```
ipython3 kernel install --user --name=pryngles-lab
```

Edit the `kernel.json` file to point to the virtual environment:

```
emacs -nw /Users/jorgezuluagacallejas/Library/Jupyter/kernels/pryngles-lab/kernel.json
```

And change the first variable of the `argv` placing there the
directory of the virtual environment:

```
{
 "argv": [
  "/Users/jorgezuluagacallejas/.virtualenvs/pryngles-lab/bin/python",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "pryngles-lab",
 "language": "python",
 "metadata": {
  "debugger": true
 }
}
```

Alternatively you can use `virtualenv`.

```
virtualenv pryngles-lab
```

Once created go to that directory and activate the virtual environment:

```
cd pringles-lab
source bin/activate
```

Install pryngles from test repo
-------------------------------

Change to a virtual environment:

```
workon pryngles-lab
```

Install dependencies:

```
python -m pip install rebound scipy ipython matplotlib tqdm dill spiceypy cmasher
```

To install pryngles from the test repo use:

```
python -m pip install --index-url https://test.pypi.org/simple/ pryngles==x.y.z.w
```

Install there a new version of pryngles (or update it).  Test if installation was succesful:

```
python -c "from pryngles import *"
```

and test the version installed:

```
python -c "from pryngles import __version__;print(__version__)"
```

To leave the virtual environmente use:

```
deactivate
```

Create a branch 
---------------

If you want to create a branch `newbranch` use:

```
git checkout -b newbranch
```

Change to another branch:

```
git checkout master
```

Fuentes
-------

- https://janakiev.com/blog/jupyter-virtual-envs/

- https://ripon-banik.medium.com/jupyter-notebook-is-unable-to-find-module-in-virtual-environment-fa0725c3f8fd

- https://virtualenvwrapper.readthedocs.io/en/latest/

- GitHub: https://git-scm.com/book/en/v2

- 
