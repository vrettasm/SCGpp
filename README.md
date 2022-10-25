# Scaled Conjugate Gradient (SCG)

Scaled Conjugate Gradient (SCG) optimization method in C++11.

### References

The original paper, that introduced this method is described in:

   1. M. F. Moller (1993). "A scaled conjugate gradient algorithm for fast
     supervised learning". Neural Networks, Volume 6, Issue 4, pp:525-533.

### Requirements

   > Boost library is required

### Examples

The examples have been compiled (successfully) on OSX10.14 with:

  > g++ --version
  >
  > Apple LLVM version 10.0.1 (clang-1001.0.46.4)
  > Target: x86_64-apple-darwin18.7.0
  >

1. [Rosenbrock](examples/example_rosenbrock.cpp), compile with:

    > g++ -std=c++11 -Wall -g ../src/common/*.cpp example_rosenbrock.cpp -o demo_rosen

2. [Sphere](examples/example_sphere.cpp), compile with:

    > g++ -std=c++11 -Wall -g ../src/common/*.cpp example_sphere.cpp -o demo_sphere

### Unittests

   Coming soon
