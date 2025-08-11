#include <Rcpp.h>
using namespace Rcpp;

#include "glm.h"


double get_sign(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

std::pair<double,bool> glm_one_group(int nlibs, const double* counts, const double* offset,
        const double* disp, const double* weights, int maxit, double tolerance, double cur_beta,
        double lambda_reg, double alpha_reg) {
    /* Setting up initial values for beta as the log of the mean of the ratio of counts to offsets.
 	 * This is the exact solution for the gamma distribution (which is the limit of the NB as
 	 * the dispersion goes to infinity. However, if cur_beta is not NA, then we assume it's good.
 	 */
	bool nonzero=false;
	if (ISNA(cur_beta)) {
		cur_beta=0;
 	   	double totweight=0;
		for (int j=0; j<nlibs; ++j) {
			const double& cur_val=counts[j];
			if (cur_val>low_value) {
				cur_beta+=cur_val/std::exp(offset[j]) * weights[j];
				nonzero=true;
			}
			totweight+=weights[j];
		}
		cur_beta=std::log(cur_beta/totweight);
	} else {
		for (int j=0; j<nlibs; ++j) {
			if (counts[j] > low_value) {
                nonzero=true;
                break;
            }
		}
	}

	// Skipping to a result for all-zero rows.
	if (!nonzero) {
        return std::make_pair(R_NegInf, true);
    }

  //double alpha = 0.5; // balance between L1 and L2 (0 = only L2, 1 = only L1)
  //double lambda = 1; // regularization strength

	// Newton-Raphson iterations to converge to mean.
  bool has_converged=false;
	for (int i=0; i<maxit; ++i) {
		double dl=0;
 	  double info=0;
		for (int j=0; j<nlibs; ++j) {
			const double mu=std::exp(cur_beta+offset[j]), denominator=1+mu*disp[j];
			dl+=(counts[j]-mu)/denominator * weights[j];
			info+=mu/denominator * weights[j];
		}

		// penalty
	  dl-=lambda_reg * (alpha_reg * get_sign(cur_beta) + (1 - alpha_reg) * cur_beta);
		info+=lambda_reg * (1 - alpha_reg);

		const double step=dl/(info);
		cur_beta+=step;
		if (std::abs(step)<tolerance) {
			has_converged=true;
			break;
		}
	}

	return std::make_pair(cur_beta, has_converged);
}


