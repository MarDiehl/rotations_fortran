# 3D Rotations

Fortran version of a 3D rotations conversion library.
Based on [marcdegraef/3Drotations](https://github.com/marcdegraef/3Drotations)

## Getting started

```
export FC=gfortran # or ifort or pgfortran
mkdir build
cd build
cmake .. -DCMAKE_Fortran_COMPILER=$FC
make
./src/test_rotations
```

## Prerequisites
- Fortran compiler
  - GNU, version 8.0 or newer (`gfortran`)
  - Intel, version 18.0 or newer (`ifort`)
  - PGI, not tested (`pgfortran`)
- cmake, version 3.10 or newer
- LAPACK library in standard location

## License

[GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Issues

- The tests for Euler angles fail, but this is because comparing Euler angles is difficult.
- Tests for some functions are missing
