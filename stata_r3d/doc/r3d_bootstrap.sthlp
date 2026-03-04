{smcl}
{* *! version 1.0.0  2025}{...}
{vieweralsosee "r3d" "help r3d"}{...}
{viewerjumpto "Syntax" "r3d_bootstrap##syntax"}{...}
{viewerjumpto "Description" "r3d_bootstrap##description"}{...}
{viewerjumpto "Options" "r3d_bootstrap##options"}{...}
{viewerjumpto "Examples" "r3d_bootstrap##examples"}{...}
{title:Title}

{phang}
{bf:r3d_bootstrap} {hline 2} Bootstrap inference for R3D estimates


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:r3d_bootstrap} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt r:eps(#)}}number of bootstrap replications; default is {cmd:reps(999)}{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt t:ests(string)}}hypothesis tests to perform; default is {cmd:tests(nullity homogeneity)}{p_end}
{synopt:{opt sub:range(numlist)}}quantile subrange for testing{p_end}
{synopt:{opt seed(#)}}set random seed{p_end}
{synopt:{opt noi:sily}}display bootstrap progress{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:r3d_bootstrap} performs bootstrap inference after {cmd:r3d} estimation.
It constructs uniform confidence bands and performs hypothesis tests on
the distributional treatment effects.

{pstd}
This command must be run after {cmd:r3d}. It uses the stored estimation
results to perform multiplier bootstrap inference.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt reps(#)} specifies the number of bootstrap replications. The default
is 999.

{phang}
{opt level(#)} sets the confidence level for confidence intervals. The default
is {cmd:level(95)}.

{phang}
{opt tests(string)} specifies which hypothesis tests to perform. Options are:
{cmd:nullity} (H₀: τ(q) = 0 for all q), {cmd:homogeneity} (H₀: τ(q) is constant),
and {cmd:gini} (H₀: equal Gini coefficients). Default is {cmd:nullity homogeneity}.

{phang}
{opt subrange(numlist)} restricts hypothesis tests to a subrange of quantiles.
Specify as {cmd:subrange(0.25 0.75)} to test only the interquartile range.

{phang}
{opt seed(#)} sets the random seed for reproducible results.

{phang}
{opt noisily} displays the bootstrap progress.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. r3d x y*, cutoff(0)}{p_end}
{phang}{cmd:. r3d_bootstrap}{p_end}

{phang}{cmd:. r3d_bootstrap, reps(1999) tests(nullity)}{p_end}

{phang}{cmd:. r3d_bootstrap, subrange(0.1 0.9) level(99)}{p_end}

{title:Stored results}

{pstd}
{cmd:r3d_bootstrap} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(reps)}}number of bootstrap replications{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(ci_lower)}}lower confidence bounds{p_end}
{synopt:{cmd:r(ci_upper)}}upper confidence bounds{p_end}
{synopt:{cmd:r(pvalues)}}p-values for hypothesis tests{p_end}
{synopt:{cmd:r(test_stats)}}test statistics{p_end}
{synopt:{cmd:r(tau_boot)}}bootstrap distribution{p_end}

{title:Author}

{pstd}David Van Dijcke, University of Michigan{p_end}
{pstd}Email: dvdijcke@umich.edu{p_end}