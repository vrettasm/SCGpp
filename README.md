# Scaled Conjugate Gradient (SCG)

Scaled Conjugate Gradient (SCG) optimization method in C++11.

### References

The original paper, that introduced this method is described in:

   1. M. F. Moller (1993). "A scaled conjugate gradient algorithm for fast
     supervised learning". Neural Networks, Volume 6, Issue 4, pp:525-533.

### Requirements

   > C++11

### Examples

The examples have been compiled (successfully) on OSX10.14 with:

  > g++ --version
  >
  > Apple LLVM version 10.0.1 (clang-1001.0.46.4)
  > Target: x86_64-apple-darwin18.7.0
  >

1. [Rosenbrock](examples/example_rosenbrock.cpp), compile with:

    `g++ -std=c++11 -Wall -g ../src/common/*.cpp example_rosenbrock.cpp -o demo_rosen`

    > [Note the true minimum is f(1.0, 1.0) = 0.0]:

      Rosenbrock example:

      SCG: Optimization started ...

      It=     0: F(x)=      3.98998 -:- Sum(|Gradients|)=      2.01935 -:- Delta(Elapsed)= 14 μsec.

      It=    50: F(x)=      1.44202 -:- Sum(|Gradients|)=      18.0427 -:- Delta(Elapsed)= 310 μsec.

      It=   100: F(x)=     0.244361 -:- Sum(|Gradients|)=     0.982217 -:- Delta(Elapsed)= 301 μsec.

      It=   150: F(x)=    0.0689755 -:- Sum(|Gradients|)=     0.380681 -:- Delta(Elapsed)= 347 μsec.

      It=   200: F(x)=    0.0445745 -:- Sum(|Gradients|)=     0.286013 -:- Delta(Elapsed)= 345 μsec.

      It=   250: F(x)=    0.0301431 -:- Sum(|Gradients|)=      0.22606 -:- Delta(Elapsed)= 384 μsec.

      It=   300: F(x)=    0.0199049 -:- Sum(|Gradients|)=     0.179038 -:- Delta(Elapsed)= 352 μsec.

      It=   350: F(x)=    0.0118292 -:- Sum(|Gradients|)=     0.135607 -:- Delta(Elapsed)= 345 μsec.

      It=   400: F(x)=   0.00479381 -:- Sum(|Gradients|)=    0.0851121 -:- Delta(Elapsed)= 345 μsec.

      It=   450: F(x)=  6.59135e-06 -:- Sum(|Gradients|)=    0.0032239 -:- Delta(Elapsed)= 347 μsec.

      It=   500: F(x)=  3.68206e-06 -:- Sum(|Gradients|)=   0.00250468 -:- Delta(Elapsed)= 346 μsec.

      It=   550: F(x)=  2.82877e-06 -:- Sum(|Gradients|)=   0.00225179 -:- Delta(Elapsed)= 345 μsec.

      Minimum 2.63299e-06 found at f(0.998378, 0.996754).

2. [Sphere](examples/example_sphere.cpp), compile with:

    `g++ -std=c++11 -Wall -g ../src/common/*.cpp example_sphere.cpp -o demo_sphere`

    > [Note the true minimum is f(0.0, 0.0, 0.0, 0.0) = 0.0]:

      Sphere example:

      SCG: Optimization started ...

      It=     0: F(x)=      410.333 -:- Sum(|Gradients|)=           70 -:- Delta(Elapsed)= 16 μsec.

      It=     2: F(x)=     0.202634 -:- Sum(|Gradients|)=      1.55556 -:- Delta(Elapsed)= 44 μsec.

      It=     4: F(x)=  6.43852e-07 -:- Sum(|Gradients|)=   0.00277283 -:- Delta(Elapsed)= 29 μsec.

      It=     6: F(x)=  9.15756e-15 -:- Sum(|Gradients|)=  3.30689e-07 -:- Delta(Elapsed)= 29 μsec.

      Minimum 9.15756e-15 found at f(1.88965e-08, 3.46436e-08, 8.18849e-08, -2.99195e-08).

### Unittests

   - [ ] Coming soon
