#include <iostream>
#include <utility>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

// Custom code.
#include "../src/include/scaled_conjugate_gradient.hpp"

class Shpere: public OptimizationSpace::OptimizationFunction {
  /*
    The Shpere function is used as a performance test problem
    for many optimization algorithms. */

  public:

    /** @brief Constructor with no parameters. */
    Shpere() {};

    /** @brief Return the function value f(x). */
    std::pair<double, std::vector<double>> operator()(std::vector<double>& x) {

      // Function value.
      double fx = 0.0;

      // Sum_i=1^n x^2_i
      for (auto& xi: x) {
        fx += std::pow(xi, 2);
      }

      // Calculate the gradient vector.
      std::vector<double> dfx(x);
      std::transform(dfx.begin(), dfx.end(), dfx.begin(),
                     std::bind(std::multiplies<double>(),
                               std::placeholders::_1, 2.0));

      // Return the f(x, y)
      return std::make_pair(fx, dfx);
    }
};


int main(int argc, char* argv[]) {

  // Display info.
  std::cout << " Sphere example: " << std::endl;

  // Create a Shpere object with n=3.
  Shpere shere_func;

  // Hamiltonian MC object.
  OptimizationSpace::SCG scaled_optim(shere_func);

  scaled_optim.set_max_it(100);
  scaled_optim.set_upd_frequency(2);

  // Initial search point.
  std::vector<double> x0{12.0, 22.0, 52.0, -19.0};

  // Start the optimization process.
  std::tuple<std::vector<double>, double> result = scaled_optim.run(x0);

  // Extract the optimal values.
  std::vector<double> x_opt = std::get<0>(result);
  double f_opt = std::get<1>(result);

  // Print the minimum values.
  std::cout << "\nMinimum " << f_opt << " found at f(";

  for (size_t i = 0; i < x_opt.size(); ++i) {

    if (i == x_opt.size()-1) {

      std::cout << x_opt[i] << ")";

    } else {

      std::cout << x_opt[i] << ", ";

    }

  }

  std::cout << std::endl;

  // Exit.
  return 0;
}
