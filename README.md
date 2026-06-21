## What is this

Botsch-Kobbelt remesher implemented using libigl and eigen, supporting arbitrary vertex attributes and selections. Run via python bindings (made with nanobind).

This is a fork of
[this repo](https://github.com/yoterel/adaptive_isotropic_remeshing), adding interpolation of arbitrary vertex-based quantities as well as support for restricting remeshing to a selection of vertices.

That repo is itself a fork of
[this repo](https://github.com/sgsellan/botsch-kobbelt-remesher-libigl), adding adaptive remeshing.

## Building

```sh
# First activate the python environment to build for; must be python>=3.8 and have numpy installed
conda activate myenv

# then build the extension
cmake -S . -B build
cmake --build build
```

A successful build will yield a shared library file in the `build` directory. On Linux, this will have a filename that looks like `bkremeshlerps.cpython*.so` (with a python version string corresponding to the python version in your environment).

In python code, you can import this shared library file as a module:
```
sys.path.insert(0, "build")
import bkremeshlerps
```

## Usage

More polished documentation is todo, but the auto-generated docs from nanobind you get by running `help(bkremeshlerps)` at a python REPL should contain all the information you need (type signatures and docstrings defined in the C++ source).

At a glance:

```python
vattrsr, fr, vselr, fiprojv, new2oldFi = bkremeshlerps.remesh_botsch_with_interps(
    vattrs,
        # (n_verts, n_features)
    f,
        # (n_faces, 3)
    vsel,
        # (n_verts,) float, where 1 means selected and remeshing allowed, 0 means not
    targetlen=avglen,
        # target edge length for remeshing. can be a single float, or a 'sizing field' (n_verts,), roughly specifying the desired surrounding edge lengths
    selection_threshold=0.5,
        # threshold for considering a vertex selected (for interpolating the selection field between iters)
    iterations=2,
        # number of iterations
    project=True,
        # whether to project final vertices back to the original surface
    smooth=True,
        # smooth should basically always be True
    verbose=True,
)

# Returns remeshed mesh, with
# vattrsr: (new_n_verts, n_features) interpolated features on remeshed vertices
# fr: (new_n_faces,) new faces
# vselr: (new_n_verts,) interpolated selection
# fiprojv: (new_n_verts,) original face indices containing vertices projected back onto the source mesh
# new2oldFi: (new_n_faces,) a best-effort mapping of the faces in fr to the original faces f. (-1) means a face was newly created or otherwise has no obvious single parent in f. (experimental, subject to change!)
```
or adaptive remeshing:
```python
vattrsr, fr, vselr, fiprojv, new2oldFi = bkremeshlerps.remesh_botsch_adaptive_with_interps(
    vattrs,
    f,
    vsel,
    epsilon=0.001,
    adaptive=True,
    selection_threshold=0.5,
    iterations=3,
    project=True,
    smooth=True,
    verbose=True,
)
```