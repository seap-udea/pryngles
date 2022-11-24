Guide to developers
===================

In this guide you will find common procedures for developers.

NOTE:
- You will need /bin/bash installed for Makefile properly works.

Run pixx on mac:
---------------

- Install intel homebrew:

  ```
  mkdir homebrew
  curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew
  mv homebrew /usr/local
  ```

- Add `/usr/local/bin` to $PATH in `.bashrc`.

- Create an alias:

  ```
  alias axbrew='arch -x86_64 /usr/local/homebrew/bin/brew'
  ```

- Uninstall llvm from the current homebrew:

  ```
  brew uninstall llvm
  ```

- Install llvm from the current homebrew:

  ```
  axbrew uninstall llvm
  ```

- Check architecture:

  ```
  lipo -info src/pryngles/pixx.cpython-39-darwin.so
  ```


Contribute
----------

For the contributor:

- Fork. Creata fork of the repository.

- Clone. Clone the repository locally.

- Contribute. Change the files you need to change.

- Stash. Create a copy of the files you have changed.
  ```
  git stash
  ```

- Pull from master.
  ```
  git pull upstream <branch>
  ```

- Pop.
  ```
  git stash pop
  ```

- Commit and push.
  ```
  git commit -am "Message"
  ```

- Pull request.

For the main developer:

- Accept pull request.

- Stash.
  ```
  git stash
  ```

- Pull.
  ```
  git pull
  ```

- Pop.
  ```
  git stash pop
  ```

- Commit and push.
  ```
  git commit -am "Message"
  ```

Convert 
-------

Convert is the task in which the Jupyter notebooks are converted into python
files.  Conversion is required for releasing the packages.

There are two types of conversion:

- **Simple conversion**: Developing Jupyter notebooks are first
  converted into Python files and then parsed to extract testing
  codes.

  In simple conversion testing code starts with "if IN_JUPYTER:"
  conditional.  All code after a cell including the "--End--" string
  is ignored in simple conversion.

- **Advanced conversion**: Also called `xconvert`.  See *JupDev*
  section below. In the first stage of this procedure a simple
  conversion is performed.  Then the resulting python module files are
  parsed searching comments of the form "#@command" instructing the
  script to separate componentes by classes and constants.

  All classes in a file are first located in a `<module>-<class>.temp`
  file.

  All stand-alone methods, which are tagged with a `#@method <class>`
  command, are then appended to the `<module>-<class>.temp` file.

  Class documentation is also identified (variables `<class>_doc`) and
  placed in a proper plass in the class file.

  Once all the classes and methods are gathered the final python
  modules are assambled.

  All comments starting with `# #` and the blank lines surrounding
  them, which are titles in the Jupyter files are deleted.  As a
  result, only the internal comments are preserved.

  Constants of the modules are all gathered in the `consts.py` module.

  In all cases the `#@end:<section>` comment are used to delimitate block of
  code.

JupDev
------

`JupDev` is a programming style designed by Jorge I. Zuluaga and based
on Jupyter. It uses notebooks to write down and test the modules of a
package and python and bash scripts to generate programatically the
production source codes.

The module files in this style are called `dev/pryngles-module.ipynb`.
They have all the structure of Jupyter notebooks, they could have
titles and text.

Special comments helps to identify which code will be included in the
final version of the production files:

- `#@external`: show the beginning of the external packages required
  for the specific modules.  All the packages will be compiled in a
  file __packages.temp.

- `#@end:<section>`: marks the end of a code selection.

  NOTE: The end of a module is marked by `#@end:module`.

- `#@consts:<Module>`: show the beginning of constants associated to a
  module __packages.temp.

- `#@method:<Class>`: marks the begining of a method for a given
  class.

- `#@standalone`: marks the begining of standalone code in a module.
  The standalone code is always at the end of the final realease code.

- `#@data`: marks the beginning of data.

- `#@doc:<Class>`: marks the beginning of the class doc.

- `Class`: marks the beginning of a class of the package.  A file with
  the name `__<module>_<Class>.temp` is used to compile the methods
  used afterwards.

Temporal files:

- `__<module>_<Class>.temp`: temporal file with a Class and all the
  methods of a module.

- `__consts.temp`: temporal file containing all constants in the
  package (uppercase variables)

- `__<module>.temp`: file with the standalone code.

Once the `.temp` files are generated the production source code is
obtained by merging the corresponding files in the
`src/pryngles/<module>.py` and the
`src/pryngles/tests/test-<module>.py` directories.

Error types
-----------

See https://www.tutorialsteacher.com/python/error-types-in-python

Developing Cycle
----------------

Tasks:

- Change an existing module.
- Create a new module.
- Release and test a test version.
- Release and test a production version.

### Create a new module

1. Create a module from template:

   ```
   cp dev/template.ipynb dev/pryngles-module.ipynb
   ```

2. Open module notebook in Jupyter.

3. Develop module class and methods.

4. Write module class tests:

   ```python
   if IN_JUPYTER:
       def test_fun(self):
            """
      	    self.assertEqual(self.P.Nr,8,True)
            self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
            self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
	    """

	class Test(unittest.TestCase):pass
        Test.test_fun=test_fun
	unittest.main(argv=['first-arg-is-ignored'],exit=False)
   ```

5. Convert module notebooks to source code:

   ```
   make convert
   ```

6. Test the module:

   ```
   make test MOD=<module>
   ```

7. Test the whole package:
   
   ```
   make test
   ```

### Change an existing module

1. Change the corresponding module notebook `dev/pryngles-<module>.ipynb`.

2. Convert module notebooks to source code:

   ```
   make convert
   ```

3. Test the module:

   ```
   make test MOD=<module>
   ```

4. Test the whole package:
   
   ```
   make test
   ```

### Release and test a test version

1. Check the latest version: `tail .versions`

2. Execute

   ```
   make release RELMODE=test VERSION=x.y.z.w
   ```
   user: seap-udea, pass: Fmunu;mu=0
 
   For version numbers see below.

3. Change to the test virtual environment:

   ```
   workon pryngles-lab
   ```

4. Install all new dependencies. Example:

   ```
   python -m pip install anytree celluloid
   ```

5. Update the local version:

   ```
   python -m pip install --index-url https://test.pypi.org/simple/ pryngles==x.y.z.w
   ```
   
   A more complete version is:

   

6. Open `Jupyter`:

   ```
   jupyter-notebook
   ```

7. Execute all notebooks:
   - `examples/pryngles-tutorial-quickstart.ipynb`
   - `papers/bright-side/pryngles-paper-figures.ipynb`

### Release and test a production version

1. Check the latest version: `tail .versions`

2. Execute

   ```
   make release RELMODE=release VERSION=x.y.z
   ```
   user: jorge.zuluaga, pass: Fmunu;mu=0

   For version numbers see below.

3. Change to the test virtual environment:

   ```
   workon pryngles-stable
   ```

4. Update the local version:

   ```
   python -m pip install -U pryngles
   ```

5. Open `Jupyter`:

   ```
   jupyter-notebook
   ```

6. Execute all notebooks:
   - `examples/pryngles-tutorial-quickstart.ipynb`
   - `papers/bright-side/pryngles-paper-figures.ipynb`


7. Update public repository:

   ```
   make public
   ```

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


Each time a new version is released, either in the production or test
platform, the src/pryngles/version.py is updated.

Version variable can be upload as:

```
from pryngles import *
print(version)
```

or the preferred:

```
import pryngles as pr
print(pr.version)
```

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
workon pryngles-lab
```

Alternatively, you can work with venv of python:

```
python3 -m venv pryngles-lab
```

Activate it:
```
source pryngles-lab/bin/activate
```

Install Jupyter in the new virtual environment:

```
python -m pip install jupyter
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

Now you can start jupyter in the new virtual environment and when
running a given Notebook you must recall to switch kernel.

To remove a virtualenvironment use:

```
rmvirtualenv pryngles-lab
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
python -c "import pryngles as pr;print(pr.version)"
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

- Virtual Environments in Python: https://janakiev.com/blog/jupyter-virtual-envs/

- Virtual Environments in Python: https://ripon-banik.medium.com/jupyter-notebook-is-unable-to-find-module-in-virtual-environment-fa0725c3f8fd

- Virtual Environments in Python: https://virtualenvwrapper.readthedocs.io/en/latest/

- GitHub: https://git-scm.com/book/en/v2

- Readthedocs: https://docs.readthedocs.io/en/stable/tutorial/
