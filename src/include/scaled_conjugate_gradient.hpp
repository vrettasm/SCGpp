#ifndef __SCALED_CONJUGATE_GRADIENT_HPP_
#define __SCALED_CONJUGATE_GRADIENT_HPP_

#include <tuple>
#include <string>
#include <vector>
#include <utility>

/** @brief Scaled Conjugate Gradient C++ optimizer.

    @author Michalis Vrettas, PhD.
            E-mail: vrettasm@gmail.com */
namespace OptimizationSpace {

  /** @brief Base class that provides the basic interface
      to use SCG for the optimzation problems.

      @author Michalis Vrettas, PhD.
              E-mail: vrettasm@gmail.com */
  class OptimizationFunction {
    public:
      // Constructor.
      OptimizationFunction() {};

      // Destructor.
      virtual ~OptimizationFunction() {};

      /** @brief Overloaded operator() calls directly the function
          to be optimized for a given set of model paramerers.
          
          @return the function and gradient values [f(x), df(x)] in
          a tuple.
          
          @param x: the state vector where the optimization algorithm
          will evaluate the [f(x), df(x)]. */
      virtual std::pair<double, std::vector<double>> operator()(std::vector<double>& x) = 0;
  };

  /** @brief Class which implements the SCG sampler.

      @cite[SCG] M. F. Moller (1993).
      "A scaled conjugate gradient algorithm for fast supervised learning".
      Neural Networks, Volume 6, Issue 4, pp:525-533.

      @author Michalis Vrettas, PhD.
              E-mail: vrettasm@gmail.com */
  class SCG {
    public:
      /** @brief Constructor for the SCG.
 
          @param func: Reference to an instance of a class
          that implements the optimization method. */
      SCG(OptimizationFunction& func);

      /** @brief Set the error tolerance in 'x'. */
      void set_x_tol(const double);

      /** @brief Set the error tolerance in 'fx'. */
      void set_f_tol(const double);

      /** @brief Set the maximum number of iterations. */
      void set_max_it(const int);

      /** @brief Set the display information flag. */
      void set_display(const bool);

      /** @brief Set the update frequency. */
      void set_upd_frequency(const int);

      /** @return Get the error tolerance in 'x'. */
      inline double get_x_tol() const { return x_tol; }

      /** @return Get the error tolerance in 'f(x)'. */
      inline double get_f_tol() const { return f_tol; }

      /** @return Get the maximum number of iterations. */
      inline int get_max_it() const { return max_it; }

      /** @return Get the display information flag. */
      inline bool get_display() const { return display; }

      /** @return Get the update frequency. */
      inline int get_upd_frequency() const { return upd_frequency; }

      /** @brief Launch SCG optimization procedure from starting position x0. */
      std::tuple<std::vector<double>, double> run(const std::vector<double>& x0);

  private:
  
      /** Error tolerance in 'x'. */
      double x_tol;

      /** Error tolerance in 'f(x)'. */
      double f_tol;

      /** Maximum number of iterations. */
      int max_it;

      /** Display information flag. */
      bool display;

      /** Update frequency for display information. */
      int upd_frequency;
      
      /** Reference for the class that computes the
          optimzation function. */
      OptimizationFunction& func;
      
      /** @note This will prevent accidental copy of the object. */
      SCG(const SCG&) = delete;

      /** @note This will prevent accidental assignment of the object. */
      SCG& operator=(SCG) = delete;
  };
}

#endif

// End-of-File.