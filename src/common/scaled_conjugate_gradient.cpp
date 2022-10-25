#include <cmath>
#include <ctime>
#include <chrono>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <utility>
#include <iostream>
#include <algorithm>
#include <functional> 
#include <boost/filesystem.hpp>

#include "../include/scaled_conjugate_gradient.hpp"

namespace OptimizationSpace {
  
  SCG::SCG(OptimizationFunction& _func): x_tol(1.0e-6), f_tol(1.0e-8),
  max_it(500), display(true), upd_frequency(50), func(_func){};

  void SCG::set_x_tol(const double new_value) {
    // Sanity check: positive threshold.
    if (new_value > 0.0) {
      x_tol = new_value;
    } else {
      
      throw std::invalid_argument(" SCG::set_x_tol:"
                                  " Error tolerance in 'x' can't be negative.");
    }
  }

  void SCG::set_f_tol(const double new_value) {
    // Sanity check: positive threshold.
    if (new_value > 0.0) {
      f_tol = new_value;
    } else {
      
      throw std::invalid_argument(" SCG::set_f_tol:"
                                  " Error tolerance in 'f(x)' can't be negative.");
    }
  }

  void SCG::set_max_it(const int new_value) {
    // Sanity check: positive iterations.
    if (new_value > 0) {
      max_it = new_value;
    } else {
      throw std::invalid_argument(" SCG::set_max_it:"
                                  " Number of maximum iterations can't be negative.");
    }
  }

  void SCG::set_display(const bool new_value) {
    display = new_value;
  }

  void SCG::set_upd_frequency(const int new_value) {
    // Sanity check: positive update frequency.
    if (new_value > 0) {
      upd_frequency = new_value;
    } else {
      
      throw std::invalid_argument(" SCG::set_upd_frequency:"
                                  " Update frequency can't be negative.");
    }
  }

  std::tuple<std::vector<double>, double> SCG::run(const std::vector<double>& x0) {
    
    // Sanity check: make sure we have data.
    if (x0.empty()) {
      throw std::invalid_argument(" SCG::run:"
                                  " Starting position 'x0' is empty.");
    }
    
    // Size of the parameter vector.
    const size_t dim_x = x0.size();
    
    // Make a copy of the input vector.
    std::vector<double> x(x0);
    
    // Initial sigma.
    const double sigma0 = 1.0e-3;
    
    // Initial function/gradients value.
    std::pair<double, std::vector<double>> res_now = func(x);
    
    // Extract the pair values.
    double f_now = res_now.first;
    std::vector<double> grad_new = res_now.second;
    
    // Initialize the grad_old vector.
    std::vector<double> grad_old(grad_new);
    
    // Increase function evaluations by one.
    // _stats["func_eval"] += 1

    // Store the current values (fx / dfx).
    double f_old = f_now;
    
    // Set the initial search direction.
    // calc{d = -grad_new}
    std::vector<double> d(grad_new);
    std::transform(d.begin(), d.end(), d.begin(),
                   std::bind(std::multiplies<double>(),
                             std::placeholders::_1, -1.0));
    
    // Forces the calculation of directional derivatives.
    bool success = true;

    // Counts the number of successes.
    int count_success = 0;
    
    // Initial scale parameter.
    double beta = 1.0;

    // Lower & Upper bounds on scale (beta).
    double beta_min = 1.0e-15;
    double beta_max = 1.0e+100;
    
    // Initialization of parameters.
    double kappa = 0.0;
    double theta = 0.0;
    double mu = 0.0;
    
    // Get the machine precision constant.
    const double eps_double = std::numeric_limits<double>::epsilon();
    
    // Check for verbosity.
    if (display) {
      std::cout << "SCG: Optimization started ..." << std::endl;
    }
    
    // Get the first time.
    auto t0 = std::chrono::steady_clock::now();
    
    // Main optimization loop.
    for (size_t j = 0; j < max_it; ++j) {
      
      // Calculate 1-st and 2-nd
      // directional derivatives.
      if (success) {
        
        // Inner-product.
        mu = std::inner_product(d.begin(), d.end(), grad_new.begin(), 0.0);
        
        if (mu >= 0.0) {
          // Copy the "new" gradient to the "d".
          std::copy(grad_new.begin(), grad_new.end(), d.begin());
  
          // Change direction.
          // calc{d = -d}
          std::transform(d.begin(), d.end(), d.begin(),
                         std::bind(std::multiplies<double>(),
                                   std::placeholders::_1, -1.0));
          // Re-estimate 'mu'.
          mu = std::inner_product(d.begin(), d.end(), grad_new.begin(), 0.0);
        }

        // Compute kappa.
        kappa = std::inner_product(d.begin(), d.end(), d.begin(), 0.0);
        
        // Check for termination.
        if (kappa < eps_double) {
          // Update the statistic.
          // _stats["nit"] = j + 1
          
          // Exit from here.
          return std::make_tuple(x, f_now);
        }
        
        // Update sigma and check the gradient
        // on a new direction.
        double sigma = sigma0 / std::sqrt(kappa);
        
        // Calculate: x_plus = x + (sigma * d).
        std::vector<double> x_plus(x);
        for (size_t k = 0; k < x_plus.size(); ++k) {
          x_plus[k] += sigma*d[k];
        }
        
        // Evaluate the function at the new point.
        std::pair<double, std::vector<double>> res_plus = func(x_plus);
        
        // Get the gradient df(x_plus).
        std::vector<double> g_plus = res_plus.second;

        // Increase function evaluations.
        //_stats["func_eval"] += 1

        // Preallocate the difference vector.
        std::vector<double> g_diff(g_plus.size(), 0.0);
        
        // Compute the difference in the gradients.
        // calc{g_plus - grad_new}
        std::transform(g_plus.begin(), g_plus.end(), grad_new.begin(),
                       g_diff.begin(), std::minus<double>());
        
        // Compute theta.
        theta = std::inner_product(d.begin(), d.end(), g_diff.begin(), 0.0)/sigma;
      }
      
      // Increase effective curvature.
      double delta = theta + (beta * kappa);
      
      if (delta <= 0.0) {
        delta = beta * kappa;
        beta = beta - (theta / kappa);
      }
      
      // Update 'alpha'.
      double alpha = -(mu / delta);
      
      // Evaluate the function at a new point.
      // calc{x_new = x + (alpha * d)}
      std::vector<double> x_new(x);
      for (size_t k = 0; k < x_new.size(); ++k) {
        x_new[k] += alpha*d[k];
      }
      
      // Evaluate fx and dfx at the new point.
      // NOTE: Because we haven't accepted yet this position as the
      // next 'x' state, we use the 'g_now' to store the gradient.
      std::pair<double, std::vector<double>> res_new = func(x_new);
      double f_new = res_new.first;
      std::vector<double> g_now = res_new.second;
      
      // Note that the gradient is computed anyway.
      //_stats["func_eval"] += 1
      
      // Calculate the new comparison ratio.
      double Delta = 2.0 * (f_new - f_old) / (alpha * mu);
      if (Delta >= 0.0) {

        // Set the flag.
        success = true;

        // Update counter.
        count_success += 1;

        // Copy the new values.
        f_now = f_new;

        // Update the new search position.
        std::copy(x_new.begin(), x_new.end(), x.begin());
          
      } else {
        
        // Cancel the flag.
        success = false;

        // Copy the old values.
        f_now = f_old;

        // Update the gradient vector.
        std::copy(grad_old.begin(), grad_old.end(), g_now.begin());
      }
      
      // Total gradient: j-th iteration.
      double total_grad = 0.0;
      
      // Sum over all the absolute values.
      for (auto& grad_i : g_now) {
        total_grad += std::abs(grad_i);
      }
      
      // Store statistics.
      //_stats["fx"][j] = f_now
      //_stats["dfx"][j] = total_grad

      // Used in verbose / display mode.
      if (display && (j%upd_frequency == 0)) {

        // Timer snapshot after 'j' iterations.
        auto tj = std::chrono::steady_clock::now();
        
        // Print the current info.
        std::cout << "It= " << std::setw(5) << j
                  << ": F(x)= " << std::setw(12) << f_now
                  << " -:- Sum(|Gradients|)= " << std::setw(12)
                  << total_grad << " -:- Delta(Elapsed)= "
                  << std::chrono::duration_cast<std::chrono::microseconds>(tj - t0).count()
                  << " μsec." << std::endl;
            
        // Assign the current time to 't0'.
        t0 = tj;
      }
      
      // Check for success.
      if (success) {
        
        // Compute the 'alpha * d'.
        std::vector<double> alpha_d(d.size(), 0.0);
        std::transform(alpha_d.begin(), alpha_d.end(), d.begin(),
                       std::bind(std::multiplies<double>(),
                       std::placeholders::_1, alpha));
                       
        // Check for termination.
        if ((*std::max_element(alpha_d.begin(), alpha_d.end()) <= x_tol) &&
            ( std::abs(f_new - f_old) <= f_tol)) {
          
          // Exit.
          return std::make_tuple(x, f_new);
          
        } else {
          
            // Update variables for the new position.
            f_old = f_new;

            // Copy the "new" gradient to the "old".
            std::copy(grad_new.begin(), grad_new.end(),
                      grad_old.begin());
            
            // Evaluate function/gradient at the new point.
            std::pair<double, std::vector<double>> res_new = func(x);
            
            // Extract f(x) and df(x).
            f_now = res_new.first;
            grad_new = res_new.second;

            // If the gradient is close to zero then exit.
            if (std::inner_product(grad_new.begin(), grad_new.end(),
                                   grad_new.begin(), 0.0) <= eps_double) {

              // Close enough to exit from here.
              return std::make_tuple(x, f_now);
              
            }
        }
      
      }
      
      // Adjust beta according to comparison ratio.
      if (Delta < 0.25) {
        beta = std::min(4.0 * beta, beta_max);
      }

      if (Delta > 0.75) {
        beta = std::max(0.5 * beta, beta_min);
      }
      
      // Update search direction using Polak-Ribiere formula
      // or re-start in direction of negative gradient after
      // 'dim_x' steps.
      if (count_success == dim_x) {
        
        // Copy the "new" gradient to the "d".
        std::copy(grad_new.begin(), grad_new.end(), d.begin());
                  
        // Change direction.
        // calc{d = -grad_new}
        std::transform(d.begin(), d.end(), d.begin(),
                       std::bind(std::multiplies<double>(),
                                 std::placeholders::_1, -1.0));
        // Reset the counter.
        count_success = 0;
        
      } else {
        
        // Check the flag.
        if (success) {
          
          // Preallocate the difference vector.
          std::vector<double> g_diff(grad_new.size(), 0.0);
        
          // Compute the difference in the gradients.
          // calc{grad_old - grad_new}
          std::transform(grad_old.begin(), grad_old.end(), grad_new.begin(),
                         g_diff.begin(), std::minus<double>());

          // Update gamma.
          double gamma = std::max(std::inner_product(grad_new.begin(), grad_new.end(),
                                                     g_diff.begin(), 0.0) / mu, 0.0);
          // Update direction.
          // calc{d = (gamma * d) - grad_new}
          for (size_t k = 0; k < d.size(); ++k) {
            d[k] = (gamma * d[k]) - grad_new[k];
          }
          
        }
        
      }
      
    }
    
    // Display a final (warning) to the user.
    std::cout << "SGC: Maximum number of iterations ("
              << max_it << ") has been reached." << std::endl;

    // Exit from here.
    return std::make_tuple(x, f_old);
  }
  
}
// End-of-File.