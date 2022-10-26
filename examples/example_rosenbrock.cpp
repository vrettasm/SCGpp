#include <iostream>
#include <utility>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

// Custom code.
#include "../src/include/scaled_conjugate_gradient.hpp"

class Rosenbrock: public OptimizationSpace::OptimizationFunction {
  /*
   * The Rosenbrock function is a non-convex function which is used
   * as a performance test problem for many optimization algorithms. */

  public:

    /** @brief Constructor with input parameters. */
    Rosenbrock(double a, double b): alpha(a), beta(b) {};

    /** @brief Return the function value f(x, y)
        for the given parameters alpha, beta. */
    std::pair<double, std::vector<double>> operator()(std::vector<double>& v) {

      // Extract the (x, y) values from the vector.
      double x = v[0];
      double y = v[1];

      // Pre-compute x^2.
      double x_square = std::pow(x, 2);

      // Function value.
      double fx = std::pow(alpha - x, 2) + beta*std::pow((y - x_square), 2);

      // Calculate the gradient vector.
      std::vector<double> dfx {2.0 * (x - alpha) - 4.0 * beta * x * (y - x_square),
                               2.0 * beta * (y - x_square)};
      // Return the f(x, y)
      return std::make_pair(fx, dfx);
    }

  private:

    /* Input parameter 'a'. */
    double alpha;

    /* Input parameter 'b'. */
    double beta;
};


int main(int argc, char* argv[]) {

  // Display info.
  std::cout << " Rosenbrock example: " << std::endl;

  // Create a rosenbrock object with default values.
  Rosenbrock rosen_func(1.0, 100.0);

  // Scaled conjugate object.
  OptimizationSpace::SCG scaled_optim(rosen_func);

  // Set up parameters.
  scaled_optim.set_max_it(5000);

  // Initial search point.
  std::vector<double> x0{-1.0, 1.0};

  // Start the optimization process.
  std::tuple<std::vector<double>, double> result = scaled_optim.run(x0);

  // Extract the optimal values.
  std::vector<double> x_opt = std::get<0>(result);
  double f_opt = std::get<1>(result);

  // Print the minimum values.
  std::cout << "\nMinimum " << f_opt << " found at f("
            << x_opt[0] << ", " << x_opt[1] << ")." << std::endl;

  // Exit.
  return 0;
}
