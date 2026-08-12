# Cluster number counts — `fast_mass`, inline selection function (C++)

**Status: reference re-expression, inline-selection variant** (added
2026-08-12). Production remains `NumCountsSel.so` (tabulated
`sel_function/S_stack`). Built as `NumCountsFastMass.so`.

A thin des_y3 driver holding `y3_cluster::SelFunctionInlineCore_t` (see
`src/models/sel_function_inline_t.hh`) instead of the production
`nosel_gl_detail::SelGLCore` — the Gauss-Legendre quadrature over true
richness lambda^tr (the integral `sel_function.py` performs in Python,
~197 ms/sample, the single largest line item in the reference pipeline)
is computed directly here, once per (lnM, z) GL cell, from the existing
immutable models `MOR_HOD_t`, `RichnessKernel_t`, `PlobLtrEMG_t`,
`richness_zkernel`. No separate `sel_function` module run is required for
this observable.

Output section `numcounts_fast_mass/vals` — same section the Python
fast_mass reference uses (the two are interchangeable, never run both in
one pipeline); production's own section (`numcountssel`) is different, so
this driver can co-run with production in the same pipeline for direct
comparison.

Validation: pending (see the des_y3 fast_mass validation plan — compare
`numcounts_fast_mass/vals` against `numcountssel/vals` and against the
adaptive `full_ltmz` reference on the real pipeline fiducial point;
numbers to be recorded here once measured).
