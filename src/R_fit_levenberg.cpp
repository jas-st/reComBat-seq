// Adapted from edgeR 4.0.1 (R_fit_levenberg.cpp)
// Modified by Zhasmina Stoyanova, 2026
// Changes: added lambda_reg, alpha_reg, num_threads parameters;
//          OpenMP parallelization with thread-local workspaces;
//          thread-safe copying of compressed_matrix rows via critical section.
// Original authors: edgeR team (see edgeR package).
#include "glm.h"
#include "objects.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// THESE MUST BE AT THE TOP LEVEL (GLOBAL SCOPE)
extern "C" {
    int openblas_get_num_threads(void) __attribute__((weak));
    void openblas_set_num_threads(int num_threads) __attribute__((weak));
}

SEXP fit_levenberg (SEXP y, SEXP offset, SEXP disp, SEXP weights, SEXP design,
                    SEXP beta, SEXP tol, SEXP maxit, SEXP lambda_reg, SEXP alpha_reg,
                    SEXP num_threads) {
    BEGIN_RCPP
     
    Rprintf("test");
    // OPENBLAS DEBUG
    if (openblas_get_num_threads) {
      Rprintf("Inside C++, OpenBLAS is using %d threads.\n", openblas_get_num_threads()); }
    
 
    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();

    double lamb=check_numeric_scalar(lambda_reg, "reg strength lambda");
    double alpha=check_numeric_scalar(alpha_reg, "reg control alpha");
    double n_threads=check_numeric_scalar(num_threads, "number of threads");

    // Getting and checking the dimensions of the arguments.
    Rcpp::NumericMatrix X=check_design_matrix(design, num_libs);
    const int num_coefs=X.ncol();

    Rcpp::NumericMatrix Beta(beta);
    if (Beta.nrow()!=num_tags || Beta.ncol()!=num_coefs) {
        throw std::runtime_error("dimensions of beta starting values are not consistent with other dimensions");
    }

    // Initializing pointers to the assorted features.
    compressed_matrix allo=check_CM_dims(offset, num_tags, num_libs, "offset", "count");
    compressed_matrix alld=check_CM_dims(disp, num_tags, num_libs, "dispersion", "count");
    compressed_matrix allw=check_CM_dims(weights, num_tags, num_libs, "weight", "count");

    // Setting up scalars.
    int max_it=check_integer_scalar(maxit, "maximum iterations");
    double tolerance=check_numeric_scalar(tol, "tolerance");

    // Initializing output objects.
    Rcpp::NumericMatrix out_beta(num_tags, num_coefs);
    Rcpp::NumericMatrix out_fitted(num_tags, num_libs);
    Rcpp::NumericVector out_dev(num_tags);
    Rcpp::IntegerVector out_iter(num_tags);
    Rcpp::LogicalVector out_conv(num_tags);

    std::vector<double> current(num_libs), tmp_beta(num_coefs), tmp_fitted(num_libs);

    //PARALLEL
    #ifdef _OPENMP
        if (n_threads > 0) {

          if (openblas_set_num_threads) {
              openblas_set_num_threads(1);
              Rprintf("Inside C++, OpenBLAS is using %d threads after the command.\n", openblas_get_num_threads());
          }

          omp_set_num_threads(n_threads);
          //std::cout << "Max threads available: " << omp_get_max_threads() << std::endl;
          Rprintf("Max threads available: %d\n", omp_get_max_threads());

    #pragma omp parallel
    {
      // Each thread gets its own workspace
      glm_levenberg thread_glbg(num_libs, num_coefs, X.begin(), max_it, tolerance, lamb, alpha);
      std::vector<double> thread_current(num_libs);
      std::vector<double> thread_tmp_beta(num_coefs);
      std::vector<double> thread_tmp_fitted(num_libs);

      // Add thread-local buffers for row data
      std::vector<double> thread_offset(num_libs);
      std::vector<double> thread_disp(num_libs);
      std::vector<double> thread_weights(num_libs); 

    #pragma omp single
    {
      //std::cout << "Threads running: " << omp_get_num_threads() << std::endl;
      Rprintf("Threads running: %d\n", omp_get_num_threads());
    }

    #pragma omp for schedule(static)
    for (int tag=0; tag<num_tags; ++tag) {
      counts.fill_row(tag, thread_current.data());

      auto beta_row = Beta.row(tag);
      std::copy(beta_row.begin(), beta_row.end(), thread_tmp_beta.begin());


      // Copy row data to thread-local buffers under protection
      // compressed_matrix::get_row() returns a pointer to a shared member variable
      // that gets overwritten on each call, so we must copy to thread-local storage
      // before passing to fit() i think
      // Pre-buffering all rows upfront
      // was rejected due to memory cost at scale (O(num_tags * num_libs))
    #pragma omp critical(get_rows)
    {
      const double* offset_ptr = allo.get_row(tag);
      std::copy(offset_ptr, offset_ptr + num_libs, thread_offset.begin());

      const double* disp_ptr = alld.get_row(tag);
      std::copy(disp_ptr, disp_ptr + num_libs, thread_disp.begin());

      const double* weights_ptr = allw.get_row(tag);           
      std::copy(weights_ptr, weights_ptr + num_libs, thread_weights.begin()); 
    }

    if (thread_glbg.fit(thread_current.data(),
                        thread_offset.data(),
                        thread_disp.data(),
                        thread_weights.data(),
                        thread_tmp_fitted.data(),
                        thread_tmp_beta.data())) {

    #pragma omp critical
    {
      std::stringstream errout;
      errout << "solution using Cholesky decomposition failed for tag " << tag+1;
      throw std::runtime_error(errout.str());
    }
    }

    std::copy(thread_tmp_fitted.begin(), thread_tmp_fitted.end(), out_fitted.row(tag).begin());
    std::copy(thread_tmp_beta.begin(), thread_tmp_beta.end(), out_beta.row(tag).begin());
    out_dev[tag] = thread_glbg.get_deviance();
    out_iter[tag] = thread_glbg.get_iterations();
    out_conv[tag] = thread_glbg.is_failure();
    }
    }
        }else
    #endif

          // NO PARALLEL
    {
      glm_levenberg glbg(num_libs, num_coefs, X.begin(), max_it, tolerance, lamb, alpha);
      std::vector<double> current(num_libs);
      std::vector<double> tmp_beta(num_coefs);
      std::vector<double> tmp_fitted(num_libs);

      for (int tag=0; tag<num_tags; ++tag) {
        counts.fill_row(tag, current.data());
        auto beta_row = Beta.row(tag);
        std::copy(beta_row.begin(), beta_row.end(), tmp_beta.begin());

        if (glbg.fit(current.data(), allo.get_row(tag), alld.get_row(tag),
                     allw.get_row(tag), tmp_fitted.data(), tmp_beta.data())) {
          std::stringstream errout;
          errout << "solution using Cholesky decomposition failed for tag " << tag+1;
          throw std::runtime_error(errout.str());
        }

        std::copy(tmp_fitted.begin(), tmp_fitted.end(), out_fitted.row(tag).begin());
        std::copy(tmp_beta.begin(), tmp_beta.end(), out_beta.row(tag).begin());
        out_dev[tag] = glbg.get_deviance();
        out_iter[tag] = glbg.get_iterations();
        out_conv[tag] = glbg.is_failure();
      }
    }

    return Rcpp::List::create(out_beta, out_fitted, out_dev, out_iter, out_conv);
    END_RCPP
}
