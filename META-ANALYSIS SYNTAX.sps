* Encoding: UTF-8.
*The below syntax creates new variables and calls on macros created by David Wilson to conduct the meta-analyses.

*Calculate standard error.
COMPUTE V=1 / (SQRT(N - 3)).
*Calculate inverse variance weight.
COMPUTE W = 1/V.
*Calculate additional variables for additional calculations.
COMPUTE WES = W*D.
COMPUTE WESSQ = W*D**2.
COMPUTE WSQ = W**2.
EXECUTE.

VARIABLE LABELS
V 'Inverse variance'
W 'Inverse variance weight'
WES 'Weighted mean effect size'
WESSQ 'Weighted mean effect size squared'
WSQ 'Inverse variance weight squared'.

DESCRIPTIVES
    /VARIABLES W WES WESSQ WSQ
    /STATISTICS SUM.

*Call the fixed/random effects macro.
INCLUDE "[macro location]\MEANES.SPS" .
MEANES ES = D / W = W /PRINT=IVZR.
MEANES ES = D / W = W.
EXECUTE.

*Call the mixed effect macro.
INCLUDE "[macro location]\METAF.SPS" .
METAF ES=D /W=W / GROUP = GENDER /MODEL=ML /PRINT.
METAF ES = D /W = W /GROUP = RACE /MODEL=ML /PRINT.
METAF ES = D /W = W /GROUP = LOCATION /MODEL=ML /PRINT.
METAF ES = D /W = W /GROUP = DESIGN /MODEL=ML /PRINT.


