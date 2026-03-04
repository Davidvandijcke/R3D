# R3D Stata Hand-Off Manual

## Goal
Finish the Stata/Mata backend so the public commands in `ado/` match the R reference implementation while defaulting to the Fortran plugin (`plugin/r3d_plugin.plugin`). Deliver production-ready code—no placeholders or shortcuts.

## Key References
- Stata front ends: `ado/r3d.ado`, `ado/r3d_bwselect.ado`, `ado/r3d_bootstrap.ado`
- Mata back end: `mata/r3d_mata.mata`
- Plugin glue: `mata/r3d_plugin_interface.mata`, `ado/r3d_plugin_load.ado`
- R specification: `R/r3d_main.R`, `R/r3d_bwselect.R`, `R/r3d_bootstrap.R`, `R/r3d_utils.R`, `R/r3d_methods.R`
- Remaining drivers: `FINAL_TEST_COMMANDS.do`, `example.do`, any surviving `test_*.do` scripts (inspect repo before use)

## Guardrails
1. Fortran plugin must be used first; fall back to Mata only when the plugin returns non-zero.
2. No mock logic or hard-coded shortcuts—everything must compute the real estimator.
3. Mirror the R package exactly for syntax, defaults, numerics, and outputs.
4. Final code must compile with `mata set matastrict on`.
5. Do not modify external repos (DiSCo, etc.).
6. Keep plugin loader/function names intact.
7. Follow the existing code style (four-space indent, lowercase keywords, minimal `//` comments).

## Work Plan

1. **Reset Mata Source**  
   Copy the clean template (`ado_test/r/r3d_mata.mata` or a known-good historical version) into `mata/r3d_mata.mata`. Remove any temporary `mata set matastrict off` line.

2. **Refactor `r3d_bandwidth_select()`**  
   Declare all temporaries at the top; replicate the three-step procedure in `R/r3d_bwselect.R` for both sharp and fuzzy cases. Use plugin-first pilot fits (`r3d_plugin_locweights`) and fall back to `r3d_locpoly` only if needed.

3. **Refactor `r3d_estimate_core()`**  
   Apply the same declaration discipline. Ensure plugin-first computation for numerator/denominator fits, re-create residuals, Frechet projections, and fuzzy adjustments exactly as the R code. Confirm stored matrices match R outputs.

4. **Refactor `r3d_bootstrap()`**  
   Move loop counters (`j`, `b`, `r`, etc.) to the top, then implement multiplier bootstrap/testing exactly as `R/r3d_bootstrap.R` (uniform bands and nullity/homogeneity/gini tests). Handle `tests()` options the same way R does.

5. **Helper Functions**  
   Review `r3d_prepare_bandwidth_matrix`, `r3d_locpoly`, `r3d_isotonic`, Gini helpers, and ensure they obey matastrict rules and mirror R logic.

6. **Compile & Verify**  
   With matastrict on, run `stata-se -b do mata/r3d_mata.mata`. Execute `FINAL_TEST_COMMANDS.do`, `example.do`, and any remaining `test_*.do` scripts. If legacy tests are missing, create lightweight checks that mimic R examples (do not keep them in the repo unless requested).

7. **Documentation & Cleanup**  
   Update `doc/r3d.sthlp` for any syntax/option changes. Confirm `install_r3d.do` still works. Delete temporary logs/scripts.

## Reminders
- Work within the current directory structure.
- Keep commits focused and well-described.
- Test after each major function rewrite.
- When behaviour is unclear, re-read the R source and replicate it precisely.
- If a needed test script was removed, craft a temporary check instead of restoring old bloat.

Following these steps will bring the Mata backend back in sync with the R package while keeping the Fortran plugin front-and-center.
