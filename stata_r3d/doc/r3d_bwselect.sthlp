{smcl}
{* *! version 1.0.0  2025}{...}
{vieweralsosee "r3d" "help r3d"}{...}
{viewerjumpto "Syntax" "r3d_bwselect##syntax"}{...}
{viewerjumpto "Description" "r3d_bwselect##description"}{...}
{viewerjumpto "Options" "r3d_bwselect##options"}{...}
{viewerjumpto "Examples" "r3d_bwselect##examples"}{...}
{title:Title}

{phang}
{bf:r3d_bwselect} {hline 2} Bandwidth selection for R3D


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:r3d_bwselect} {it:runvar} {it:outcomevar1} [{it:outcomevar2} ...] {ifin}
{cmd:,} {opt c:utoff(#)} [{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt c:utoff(#)}}specifies the RD cutoff value{p_end}
{synopt:{opt m:ethod(string)}}bandwidth selection method: {cmd:simple} (default) or {cmd:frechet}{p_end}
{synopt:{opt fuz:zy(varname)}}treatment variable for fuzzy RD{p_end}
{synopt:{opt pol:ynomial(#)}}polynomial order; default is 1{p_end}
{synopt:{opt k:ernel(string)}}kernel function: {cmd:triangular} (default), {cmd:epanechnikov}, or {cmd:uniform}{p_end}
{synopt:{opt q:uantiles(numlist)}}specific quantiles to evaluate{p_end}
{synopt:{opt nq:uantiles(#)}}number of equally-spaced quantiles; default is 100{p_end}
{synopt:{opt bwg:rid(numlist)}}bandwidth grid for search{p_end}
{synopt:{opt w:eights(varname)}}weights for observations{p_end}
{synopt:{opt ver:bose}}display detailed output{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:r3d_bwselect} performs data-driven bandwidth selection for R3D estimation.
It implements MSE-optimal and IMSE-optimal bandwidth selection procedures.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt cutoff(#)} specifies the RD cutoff value. This option is required.

{phang}
{opt method(string)} specifies the bandwidth selection method. {cmd:simple} uses
MSE-optimal bandwidths averaged across quantiles. {cmd:frechet} uses IMSE-optimal
aggregation. Default is {cmd:simple}.

{phang}
{opt fuzzy(varname)} specifies the treatment variable for fuzzy RD designs.

{phang}
{opt polynomial(#)} specifies the order of the local polynomial. Default is 1.

{phang}
{opt kernel(string)} specifies the kernel function. Options are {cmd:triangular}
(default), {cmd:epanechnikov}, or {cmd:uniform}.

{phang}
{opt quantiles(numlist)} specifies particular quantiles at which to evaluate
bandwidth selection.

{phang}
{opt nquantiles(#)} specifies the number of equally-spaced quantiles. Default is 100.

{phang}
{opt bwgrid(numlist)} specifies the bandwidth grid for the search procedure.
If not specified, a grid around Silverman's rule of thumb is used.

{phang}
{opt weights(varname)} specifies weights for observations.

{phang}
{opt verbose} displays detailed output including quantile-specific bandwidths.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. r3d_bwselect score y1 y2 y3, cutoff(0)}{p_end}

{phang}{cmd:. r3d_bwselect x y*, cutoff(70) fuzzy(treatment) verbose}{p_end}

{phang}{cmd:. r3d_bwselect running y*, c(0) method(frechet) quantiles(0.1 0.5 0.9)}{p_end}

{title:Stored results}

{pstd}
{cmd:r3d_bwselect} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(h_star)}}selected bandwidth{p_end}
{synopt:{cmd:r(h_den)}}bandwidth for treatment (fuzzy RD only){p_end}
{synopt:{cmd:r(h0)}}initial bandwidth (Silverman's rule){p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(polynomial)}}polynomial order{p_end}
{synopt:{cmd:r(nquantiles)}}number of quantiles{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(method)}}bandwidth selection method{p_end}
{synopt:{cmd:r(kernel)}}kernel function{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(h_quantiles)}}quantile-specific optimal bandwidths{p_end}

{title:Author}

{pstd}David Van Dijcke, University of Michigan{p_end}
{pstd}Email: dvdijcke@umich.edu{p_end}