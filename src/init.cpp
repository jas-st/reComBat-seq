#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(compute_nbdev, 5),
	CALLDEF(compute_apl, 6),
	CALLDEF(maximize_interpolant, 2),

	CALLDEF(fit_levenberg, 10),
	CALLDEF(get_levenberg_start, 6),
	CALLDEF(fit_one_group, 9),
	CALLDEF(get_one_way_fitted, 3),
	CALLDEF(add_prior_count, 3),
    CALLDEF(ave_log_cpm, 9),

	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {NULL, NULL, 0}
  };

void attribute_visible R_init_recombatseqv2(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
