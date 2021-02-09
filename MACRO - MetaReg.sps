*--------------------------------------------------------------
*' SPSS/Win 6.1 or Higher Macro -- Written by David B. Wilson
*' Meta-Analysis Modified Weighted Multiple Regression for
*' any type of effect size
*' To use, initialize macro with the include statement:
*' INCLUDE "[drive][path]METAREG.SPS" .
*' Syntax for macro:
*' METAREG ES=varname /W=varname /IVS=varlist 
*'   /MODEL=option /PRINT=option .
*' Where ES is the effect size variable, W is the inverse
*' variance weight, IVS is the list of independent variables
*' and MODEL is either FE for a fixed effects model, MM for
*' a random effects model estimated via the method of moments,
*' and ML is a random effects model estimated via iterative
*' maximum likelihood.  If /MODEL is omitted, FE is the
*' default.  The /PRINT subcommand has the option EXP and
*' if specified will print the exponent of the B coefficient
*' (the odds-ratio) rather than beta.  If /PRINT is omitted, 
*' beta is printed.
*' Example:
*'
*' metareg es = effct /w = invweght /ivs = txvar1 txvar2
*'    /model = fe .
*'
*' Version 2005.05.23
*'
*--------------------------------------------------------------
preserve
set printback=off
define metareg (es=!charend('/')
 /w=!charend('/') /ivs=!charend('/')
 /model = !default('FE') !charend('/')
 /print = !default('RAW') !charend('/'))
preserve
set printback=off mprint off 

*--------------------------------------------------------------
* Enter matrix mode
*--------------------------------------------------------------
matrix 

*--------------------------------------------------------------
*  Get data from active file
*--------------------------------------------------------------
get data /file * /variables = !es !w !ivs /missing omit

*--------------------------------------------------------------
*  Create vectors and matrices
*--------------------------------------------------------------.
compute es = data(1:nrow(data),1) .
compute w = data(1:nrow(data),2) .
compute x = make(nrow(data),ncol(data)-1,1) .
do if ncol(data)>3 .
+ compute x(1:nrow(data),2:(ncol(data)-1))=
                  data(1:nrow(data),3:(ncol(data))).
+ else if ncol(data)=3 .
+ compute x(1:nrow(data),2)=data(1:nrow(data),3).
end if .
compute p = make(1,ncol(x),1) .
compute k = make(1,nrow(x),1) .
compute v = (w&**-1) .
release data .

*--------------------------------------------------------------
*  Recompute weights for random effects models
*  Method of moments
*--------------------------------------------------------------.
!IF (!model !eq 'MM'|!model !eq 'mm'|!model !eq 'ML'|
   !model !eq 'ml'|!model !eq 'REML'|!model !eq 'reml') !THEN .
compute xwx = T(x&*(w*p))*x .
compute B   = inv(xwx)*T(x&*(w*p))*es .
compute qe  = csum(es&*w&*es) - T(B)*xwx*B .
compute c = (qe-(nrow(es)-ncol(x)))/
     csum(w - rsum((x&*(w*p)*inv(T(x&*(w*p))*x))&*x&*(w*p))) .
do if c<0 .
+ compute c = 0 .
end if .
compute w   = 1/(v+c) .
!IFEND .

!IF (!model !eq 'ML' | !model !eq 'ml'|!model !eq 'REML'|
   !model !eq 'reml') !THEN .
compute c2 = c .
loop l=1 to 200 .
.compute loops = l .
.compute c = c2 .
.compute w = 1/(v + c) .
.compute xw = x&*(w*p) .
.compute xwx = T(xw)*x .
.compute B = inv(xwx)*T(xw)*es .
.compute r = es -x*B .
.compute c2 = csum(w&**2&*(r&**2 - v))/csum(w&**2) .
.do if c2<0 .
. compute c2 = 0 .
.end if .
end loop if abs(c2 - c)<.0000000001 .
compute c = c2 .
compute w = 1/(v + c) .
compute se_c = sqrt(2/csum(w&**2)) .
!IFEND .

!IF (!model !eq 'REML' | !model !eq 'reml') !THEN .
compute c = c2*(nrow(es)*inv(nrow(es)-ncol(x))) .
compute w = 1/(v + c) .
compute se_c = sqrt(2/csum(w&**2)) .
!IFEND .

*--------------------------------------------------------------
*  Compute Final Model
*--------------------------------------------------------------.
compute xw = x&*(w*p) .
compute xwx = T(xw)*x .
compute B = inv(xwx)*T(xw)*es .

*--------------------------------------------------------------
*  Compute Homogeneity Q for each B, for the fit of the
*  regression and for the residuals
*--------------------------------------------------------------.
compute meanes = csum(es&*w)/csum(w) .
compute q =  csum((meanes - es)&*(meanes - es)&*w) .
compute qe = csum(es&*w&*es) - T(B)*xwx*B .
compute qr = q - qe .
compute dfe = nrow(es)-ncol(x) .
compute dfr = ncol(x)-1 .
compute dft = nrow(es)-1 .
compute pe = 1 - chicdf(qe,dfe) .
compute pr = 1 - chicdf(qr,dfr) .
compute se = sqrt(diag(inv(xwx))) .
compute zvalues = B&/se .
compute pvalues = (1 - cdfnorm(abs(zvalues)))*2 .
compute lowerB = B - se*1.96 .
compute upperB = B + se*1.96 .

*--------------------------------------------------------------
*  Compute standardized coefficients (betas)
*--------------------------------------------------------------
compute d = x - T(T(csum(x&*(w*p)))*k)&/csum(w) .
compute sx = sqrt(diag(T(d&*(w*p))*d&/csum(w))) .
compute sy = sqrt(q/csum(w)) .
compute beta = (B&*sx)&/sy .
compute r2 = qr/(qr+qe) .
!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
compute beta = exp(B) .
!IFEND .

*--------------------------------------------------------------
*  Create results matrices
*--------------------------------------------------------------.
compute homog = make(3,3,-999) .
compute homog(1,1) = qr .
compute homog(1,2) = dfr .
compute homog(1,3) = pr .
compute homog(2,1) = qe .
compute homog(2,2) = dfe .
compute homog(2,3) = pe .
compute homog(3,1) = q .
compute homog(3,2) = dft .
compute homog(3,3) = 1 - chicdf(q,dft) .

compute keep = make(ncol(x),7,-999) .
compute keep(1:ncol(x),1) = B .
compute keep(1:ncol(x),2) = se .
compute keep(1:ncol(x),3) = lowerB .
compute keep(1:ncol(x),4) = upperB .
compute keep(1:ncol(x),5) = zvalues .
compute keep(1:ncol(x),6) = pvalues .
compute keep(1:ncol(x),7) = beta .

compute descrpt = make(1,3,-999) .
compute descrpt(1,1) = meanes .
compute descrpt(1,2) = r2 .
compute descrpt(1,3) = nrow(es) .

*--------------------------------------------------------------
*  Print Results
*--------------------------------------------------------------.
print /title "Version 2005.05.23".
print /title " *****  Inverse Variance Weighted Regression  *****" .
!IF (!model !eq 'ML' | !model !eq 'ml'|!model !eq 'MM'|
    !model !eq 'mm' | !model !eq 'REML' | !model !eq 'reml') !THEN .
print  /title " *****  Random Intercept, Fixed Slopes Model  *****" .
!ELSE
print  /title " *****  Fixed Effects Model via OLS  *****" .
!IFEND .

print descrpt /title "------- Descriptives -------"
  /clabel "Mean ES" "R-Squared" "k"
  /format f12.4 .
print homog /title "------- Homogeneity Analysis -------"
  /clabel "Q" "df" "p"
  /rlabel "Model" "Residual" "Total"
  /format f12.4 .
!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
print keep /title "------- Regression Coefficients -------"
    /clabel "B"  "SE" "-95% CI" "+95% CI" "Z" "P" "EXP(B)"
    /rlabel "Constant" !ivs
    /format f8.4 .
!ELSE
print keep /title "------- Regression Coefficients -------"
    /clabel "B"  "SE" "-95% CI" "+95% CI" "Z" "P" "Beta"
    /rlabel "Constant" !ivs
    /format f8.4 .
!IFEND .

!IF (!model !eq 'MM'|!model !eq 'mm') !THEN .
print c
  /title  "------- Method of Moments Random Effects Variance"+
  " Component ------- "
  /rlabel "v      ="
  /format f8.5 .
!IFEND .
!IF (!model !eq 'ML' | !model !eq 'ml') !THEN .
+ compute mlc      = make(2,1,-999) .
+ compute mlc(1,1) = c .
+ compute mlc(2,1) = se_c .
+ print mlc
   /title   "------- Maximum Likelihood Random Effects"+
   " Variance Component ------- "
   /rlabel "v      =" "se(v)  =" 
   /format f8.5 .
!IFEND .
!IF (!model !eq 'REML' | !model !eq 'reml') !THEN .
+ compute mlc      = make(2,1,-999) .
+ compute mlc(1,1) = c .
+ compute mlc(2,1) = se_c .
+ print mlc
   /title   "------- Restricted Maximum Likelihood Random"+
   " Effects Variance Component ------- "
   /rlabel "v      =" "se(v)  =" 
   /format f8.5 .
!IFEND .

*--------------------------------------------------------------
* End matrix mode
*--------------------------------------------------------------
end matrix 

*--------------------------------------------------------------
* Restore settings and exit
*--------------------------------------------------------------
restore
!enddefine
restore
