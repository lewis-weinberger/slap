# Setup

## Creating a virtualenv
We can use a virtualenv to install the dependencies in an isolated system. The following assumes you have installed Python 3 (latest version can be downloaded here: [https://www.python.org/downloads/](https://www.python.org/downloads/)), which includes the package installer [pip](https://pypi.org/project/pip/). You will also need a set of compilers for C, C++ and Fortran ([gcc](https://gcc.gnu.org/) is recommended) to compile some of the libraries.

First we install virtualenv:

```
$ pip install [--user] virtualenv
```
Note if you have an installation of Python 2 alongside your Python 3, you may need to use `pip3`. The optional `--user` flag can be used if you don't have administrative permissions for your system. We can then create our virtualenv directory (which will store the executables and libraries) using the command:

```
$ virtualenv slap-env
```
This should create a new directory called `slap-env`. To activate the virtualenv you can use the command:

```
$ source slap-env/bin/activate
```
and if successful you should see a prompt like:

```
(slap-env) $
```
We can now install the required depencies:
```
(slap-env) $ pip install --upgrade pip wheel
(slap-env) $ pip install numpy scipy
(slap-env) $ pip install hmf
(slap-env) $ pip install h5py
(slap-env) $ pip install matplotlib
```
Note that the `hmf` library uses `camb`, which needs a Fortran compiler to install. To exit the virtualenv you can use the command `deactivate`.

From within this virtualenv you should now be set up to run **slap**, an example of which is explained in [EXAMPLE](EXAMPLE.md).