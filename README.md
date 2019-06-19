# FMM_3D_pybind11

A repository providing an example of pybind11 usage with the 3D Black-box Fast Multipole Method (https://github.com/ruoxi-wang/PBBFMM3D).
The black-box FMM originates from W.Fong & E.Darve (2009).

## Requirements

- [pybind11](https://github.com/pybind/pybind11) : You need to have the path/to/pybind11/include in CPATH, or to include it manually
- [cppimport](https://github.com/tbenthompson/cppimport) : I use it, but it is not the only solution to compile the wrappers into standard (C)Python extension modules

## Usage

The consistency of the binding is evaluated by comparing the output (QH product, H a vector of unitary charges here) of the BBFMM calculated "directly" in C++ on the one hand, and in Python with pybind11 on the other hand. The former is calculated with :

```
make test_PBBFMM
cd exec
./test_PBBFMM
```
and the latter with `python3 tests_PBBFMM3D/test_PBBFMM.py`, that also diplays the results.

## Files description

### data
radialIO_3434.txt is a presaved distributions of particles (numpy arrays saved as txt). QH stands for the product of the kernel matrix Q and the sets of charges H ; these notations are often found in N-body problems literature.

### test_PBBFMM3D

  - 2Dexamples folder: contain examples I coded to help me understand pybind11 basics (pybind11::array ...)
  - PBBFMM_shot_test : source script to calculate the QH product directly with the PBBFMM3D submodule and a given distribution of particles
  - PBBFMM_binding_test.cpp and test_PBBFMM.py : pybind11 wrapper and python tests script, the distribution of particles being an array from an other Python application
