{smcl}
{* *! version 1.0.0  2025}{...}
{vieweralsosee "[R] rd" "help rd"}{...}
{vieweralsosee "[R] rdrobust" "help rdrobust"}{...}
{viewerjumpto "Syntax" "r3d##syntax"}{...}
{viewerjumpto "Description" "r3d##description"}{...}
{viewerjumpto "Options" "r3d##options"}{...}
{viewerjumpto "Examples" "r3d##examples"}{...}
{viewerjumpto "Stored results" "r3d##results"}{...}
{viewerjumpto "References" "r3d##references"}{...}
{title:Title}

{phang}
{bf:r3d} {hline 2} Regression discontinuity with distribution-valued outcomes


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:r3d} {it:runvar} {it:outcomevar1} [{it:outcomevar2} ...] {ifin}
{cmd:,} {opt c:utoff(#)} [{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt c:utoff(#)}}specifies the RD cutoff value{p_end}
{synopt:{opt m:ethod(string)}}estimation method: {cmd:simple} (default) or {cmd:frechet}{p_end}
{synopt:{opt fuz:zy(varname)}}treatment variable for fuzzy RD{p_end}
{synopt:{opt pol:ynomial(#)}}polynomial order for local regression (default is 1){p_end}
{synopt:{opt b:andwidth(numlist)}}bandwidth(s) for local regression{p_end}
{synopt:{opt bws:elect}}perform data-driven bandwidth selection{p_end}
{synopt:{opt k:ernel(string)}}kernel function: {cmd:triangular} (default), {cmd:epanechnikov}, or {cmd:uniform}{p_end}

{syntab:Quantiles}
{synopt:{opt q:uantiles(numlist)}}specific quantiles to evaluate (between 0 and 1){p_end}
{synopt:{opt nq:uantiles(#)}}number of equally-spaced quantiles (default is 100){p_end}

{syntab:Inference}
{synopt:{opt boot:strap(#)}}number of bootstrap replications for inference{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}

{syntab:Additional}
{synopt:{opt w:eights(varname)}}weights for observations{p_end}
{synopt:{opt nog:raph}}suppress graphical output{p_end}
{synopt:{opt sav:ing(filename)}}save results to dataset{p_end}
{synopt:{opt rep:lace}}overwrite existing dataset{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:r3d} implements regression discontinuity (RD) designs where the outcome
is a distribution rather than a scalar. The command estimates how the entire
distribution of the outcome variable changes at the RD cutoff.

{pstd}
The running variable is specified first, followed by one or more outcome
variables. When multiple outcome variables are specified, they are treated
as repeated measurements from the same unit's outcome distribution.

{pstd}
Two estimation methods are available:

{phang}
{cmd:simple}: Estimates treatment effects pointwise at each quantile using
local polynomial regression.

{phang}
{cmd:frechet}: Uses global Fréchet regression with isotonic projection to
ensure monotonicity of the estimated quantile function.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt cutoff(#)} specifies the RD cutoff value. This option is required.

{phang}
{opt method(string)} specifies the estimation method. {cmd:simple} (the default)
estimates treatment effects separately for each quantile. {cmd:frechet} uses
a global approach with isotonic constraints.

{phang}
{opt fuzzy(varname)} specifies the treatment variable for fuzzy RD designs.
If not specified, sharp RD is assumed.

{phang}
{opt polynomial(#)} specifies the order of the local polynomial. The default
is 1 (local linear regression).

{phang}
{opt bandwidth(numlist)} specifies the bandwidth(s) for local regression.
For sharp RD, one bandwidth is required. For fuzzy RD, two bandwidths can
be specified (for numerator and denominator).

{phang}
{opt bwselect} requests data-driven bandwidth selection using MSE-optimal
or IMSE-optimal criteria.

{phang}
{opt kernel(string)} specifies the kernel function. Options are {cmd:triangular}
(default), {cmd:epanechnikov}, or {cmd:uniform}.

{dlgtab:Quantiles}

{phang}
{opt quantiles(numlist)} specifies particular quantiles at which to evaluate
treatment effects. Values must be between 0 and 1.

{phang}
{opt nquantiles(#)} specifies the number of equally-spaced quantiles to use.
The default is 100. This option is ignored if {opt quantiles()} is specified.

{dlgtab:Inference}

{phang}
{opt bootstrap(#)} specifies the number of bootstrap replications for
constructing uniform confidence bands and hypothesis tests. If not specified,
pointwise standard errors are computed.

{phang}
{opt level(#)} sets the confidence level for confidence intervals. The default
is {cmd:level(95)}.

{dlgtab:Additional}

{phang}
{opt weights(varname)} specifies weights for observations.

{phang}
{opt nograph} suppresses the default plot of treatment effects by quantile.

{phang}
{opt saving(filename)} saves the results to a Stata dataset.

{phang}
{opt replace} overwrites an existing dataset when using {opt saving()}.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. r3d score y1 y2 y3 y4 y5, cutoff(0)}

{phang}{cmd:. r3d score y*, cutoff(70) method(frechet) bandwidth(10)}

{phang}{cmd:. r3d running_var outcome*, cutoff(0) fuzzy(treatment) bwselect}

{phang}{cmd:. r3d x y1-y10, c(0) bootstrap(999) quantiles(0.1 0.25 0.5 0.75 0.9)}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:r3d} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(cutoff)}}RD cutoff value{p_end}
{synopt:{cmd:e(polynomial)}}polynomial order{p_end}
{synopt:{cmd:e(h_num)}}bandwidth for outcome{p_end}
{synopt:{cmd:e(h_den)}}bandwidth for treatment (fuzzy RD only){p_end}
{synopt:{cmd:e(nquantiles)}}number of quantiles{p_end}
{synopt:{cmd:e(bootstrap)}}number of bootstrap replications{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:r3d}{p_end}
{synopt:{cmd:e(method)}}estimation method{p_end}
{synopt:{cmd:e(kernel)}}kernel function{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (treatment effects by quantile){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(quantiles)}}quantile grid{p_end}
{synopt:{cmd:e(pvalues)}}p-values for hypothesis tests (if bootstrap){p_end}

{marker references}{...}
{title:References}

{phang}
Van Dijcke, D. 2025. Regression Discontinuity Designs with Distributional Outcomes.
{it:Working Paper}.

{phang}
Cattaneo, M. D., Idrobo, N., and Titiunik, R. 2019. 
{it:A Practical Introduction to Regression Discontinuity Designs: Foundations}.
Cambridge University Press.

{marker author}{...}
{title:Author}

{pstd}
This Stata implementation is based on the R package by David Van Dijcke.

{pstd}
For questions and support, visit: {browse "https://github.com/dvandijcke/r3d"}