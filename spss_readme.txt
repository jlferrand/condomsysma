SPSS for Windows Meta-Analysis Macros
Written by David B. Wilson
dwilsonb@gmu.edu

Included are three macros:

	MeanES.sps
	MetaF.sps
	MetaReg.sps

To make use of these macros, you must perform your analyses from a
syntax window, rather then through the use of the pull-down menus
system.  It is also useful to change the default setting turning on
the "Display of commands in the log" found under the "Viewer" tab of
the options dialog box.

To use these macros, you must first initialize each macro with the
include statement.  This only needs to be done once during an SPSS
session.  For example, to initialize the macro MeanES.sps, run the
following SPSS syntax:

	INCLUDE 'A:\MEANES.SPS' .

The above syntax assumes that the macro is currently stored on a
floppy disk inserted into the A drive.  If it is in a different
location, modify the drive letter and path, e.g.,

	INCLUDE 'C:\PROGRAM FILES\SPSS\MACROS\MEANES.SPS' .

You can place these macros in any location that suits you.

Once initialized, a macro can be used like any other SPSS procedure.
Unfortunately, these macros are less stable than built-in features and
returns a long list of error messages if you make an error, such as
mistype a variable name.

These macros analyze the effect size as is.  That is, it is assumed
that you have performed any necessary transformations to ready the
effect size for analysis.  Common transformations are the small sample
size bias correction for d-type effect sizes, the Zr transform for
correlation (r) type effect sizes, and the natural log for odds-ratio
and risk-ratio type effect sizes.  The effect size variable specified
in the syntax statements below should be the transformed effect size.
See the section on the /PRINT statement below for information on how
to get the results printed in the original metric.

MEANES.SPS

MeanES.sps performs basic central tendency statistics, such as mean
effect size, confidence intervals, z-test, and homogeneity test.  It
calculates both fixed and random effects models.

To use this macro you need to specify the effect size variable and the
inverse variance weight variable.  The macro takes the effect size as
is, therefore, any desired transformations, e.g., Hedges' small sample
size bias correction, should already be applied.  If the adjusted
effect size variable name is "esadj" and the inverse variance weight
is "wght" then the macro syntax is:

	MEANES ES = ESADJ /W = WGHT .

METAF.SPS

MetaF.sps performs the analog to the one-way ANOVA analysis and allows
you to specify either a fixed effects model, or three methods of
estimating mixed effects (random effects) models: method-of-moments,
full-information maximum likelihood, and restricted-information
maximum likelihood.  MetaF.sps is useful for testing the differences
across mean effect sizes for a categorical variable, such as treatment
type.  Examples:

Fixed effect:
	METAF ES = ESADJ /W = WGHT /GROUP = TXTYPE .
Method of moments random effects:
	METAF ES = ESADJ /W = WGHT /GROUP = TXTYPE /MODEL = MM.
Full-information maximum likelihood:
	METAF ES = ESADJ /W = WGHT /GROUP = TXTYPE /MODEL = ML.
Restricted-information maximum likelihood:
	METAF ES = ESADJ /W = WGHT /GROUP = TXTYPE /MODEL = REML.

METAREG.SPS

MetaReg.sps performs weighted generalized least squares regression and
allows you to specify either a fixed effects model, or three methods
of estimating mixed effects (random effects) models: method-of-
moments, full-information maximum likelihood, and restricted-
information maximum likelihood.  Examples:

Fixed effect:
	METAF ES = ESADJ /W = WGHT /IVS = RANDOM TX1 TX2 .
Method of moments random effects:
	METAF ES = ESADJ /W = WGHT /IVS = RANDOM TX1 TX2 /MODEL = MM.
Full-information maximum likelihood:
	METAF ES = ESADJ /W = WGHT /IVS = RANDOM TX1 TX2 /MODEL = ML.
Restricted-information maximum likelihood:
	METAF ES = ESADJ /W = WGHT /IVS = RANDOM TX1 TX2 /MODEL = REML.

PRINT OPTIONS

MetaF.sps and MeanES.sps have an optional print statement (/PRINT) that is
useful if you are analyzing Zr transformed correlations or logged odds-
ratios.  The "/PRINT IVZR" statement converts the result back into
correlation coefficients (inverse Zr transformation).  The "/PRINT EXP"
statement prints the exponent of the results, that is, it prints odds-ratios
if the input effect sizes were logged odds-ratios.


!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!

RUNNING ALL 3 MACROS

SPSS seems unable to handle all three macros as part of a single
session.  Therefore, if you need to use all three macros you will need
to breakout your analysis run into multiple sessions (i.e., close SPSS
after you have completed using two of the macros, restart and
initialize the third macro and perform analyses).  You you find a
solution to the problem, please let me know.


