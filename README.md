utils: Utilities for VTK, VMTK and matplotlib
============================================

Author: Arjan Geers (ajgeers@gmail.com)

About
-----

This is a set of Python modules with functions I commonly use when dealing with
objects from the Visualization Toolkit ([VTK]) and the Vascular Modeling Toolkit ([VMTK]) and for making plots with [matplotlib].

Modules
-------

- `colorlib`
Module with color definitions from colorbrewer2.org
- `iolib` Module for reading, writing, compressing and converting files
- `linalglib` Module for linear algebra operations
- `plotlib` Module for plots
- `vmtklib` Module for functions in the Vascular Modeling Toolkit (VMTK)
- `vtklib` Module for operating on VTK data objects

Usage
-----

To be able to use the functions, clone the repository to your local machine and add the `utils` folder to the environment variable `PYTHONPATH`. On Mac OSX and assuming `/Users/foo/sandbox` as working directory, this can be done with:

```sh
mkdir /Users/foo/sandbox
cd /Users/foo/sandbox
git clone https://github.com/ajgeers/utils.git
export PYTHONPATH=/Users/foo/sandbox/utils/utils:$PYTHONPATH
```

To permanently add `utils` to `PYTHONPATH`, add the last line to the .bash_profile file in your home directory.

Now, in Python scripts or in a Python interpreter the functions can be called as in the following example:

```python
import iolib
import vtklib

# Slice an unstructured grid with the xy-plane at z=0,
# triangulate the slice and write it to disk.
cfd = iolib.readvtu('cfd.vtu')
xyslice = vtklib.slicedataset(cfd, [0, 0, 0], [0, 0, 1])
xyslice = vtklib.triangulate(xyslice)
iolib.writevtp(xyslice, 'xyslice.vtp')
```

Dependencies
------------

The scripts in this repository were successfully run with:
- [Python] 2.7
- [NumPy] 1.8
- [matplotlib] 1.3
- [VMTK] 1.2
- [VTK] 5.10

[Python]:http://www.python.org
[NumPy]:http://www.numpy.org
[matplotlib]:http://matplotlib.org
[VMTK]:http://www.vmtk.org
[VTK]:http://www.vtk.org
[www.vtk.org]:http://www.vtk.org
[Anaconda]:http://store.continuum.io/cshop/anaconda


License
-------

BSD 2-Clause