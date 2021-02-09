*--------------------------------------------------------------
*' SPSS/Win 6.1 or Higher Macro -- Written by David B. Wilson
*' Meta-Analysis Analog to the Oneway ANOVA for any type of ES
*' To use, initialize macro with the include statement:
*' INCLUDE "[drive][path]METAF.SPS" .
*' Syntax for macro:
*' METAF ES=varname /W=varname /GROUP=varname /MODEL=option .
*' Where ES is the effect size, W is the inverse variance 
*' weight, GROUP is the numeric categorical independent variable
*' and MODEL is either FE for a fixed effects model, MM for
*' a random effects model estimated via the method of moments,
*' and ML is a random effects model estimated via iterative
*' maximum likelihood.  If "/MODEL" is omitted, FE is the
*' default.
*' /PRINT has the options "EXP" and "IVZR".  The former
*' prints the exponent of the results (odds-ratios) and
*' the latter prints the inverse Zr transform of the 
*' results.  If the /PRINT statement is ommitted, the
*' results are printed in their raw form.

*' example:
*'
*' metaf es = effct /w = invweght /group = txvar1
*'    /model = fe .
*'
*' Version 2005.05.23
*'
*--------------------------------------------------------------
preserve
set printback=off
define metaf (es=!charend('/')
 /w=!charend('/') /group=!charend('/')
 /model = !default('FE') !charend('/')
 /print = !default('RAW') !charend('/'))
preserve
set printback=off mprint off

*--------------------------------------------------------------
* Enter matrix mode
*--------------------------------------------------------------
sort cases by !group
matrix 

*--------------------------------------------------------------
*  Get data from active file
*--------------------------------------------------------------
get data /file * /variables = !es !w !group /missing omit 

*--------------------------------------------------------------
*  Create vectors and matrices
*--------------------------------------------------------------.
compute es = data(1:nrow(data),1).
compute w = data(1:nrow(data),2).
compute v = (w&**-1) .
compute grp = data(1:nrow(data),3).
compute x = design(grp) .
compute p = make(1,ncol(x),1) .
compute k = make(1,nrow(x),1) .
compute group = inv(T(x)*x)*T(x)*grp .
release data .

*--------------------------------------------------------------
*  Recompute weights for random effects models
*  Method of moments
*--------------------------------------------------------------.
!IF (!model !eq 'MM'|!model !eq 'mm'|!model !eq 'ML'|
   !model !eq 'ml'|!model !eq 'REML'|!model !eq 'reml') !THEN .
compute xwx = T(x&*(w*p))*x .
compute B   = inv(xwx)*T(x&*(w*p))*es .
compute qw  = csum(es&*w&*es) - T(B)*xwx*B .
compute c = (qw-(nrow(es)-ncol(x)))/
     csum(w - rsum((x&*(w*p)*inv(T(x&*(w*p))*x))&*x&*(w*p))) .
do if c<0 .
+ compute c = 0 .
end if .
compute w   = 1/(v+c) .
!IFEND .

!IF (!model !eq 'ML' | !model !eq 'ml'|!model !eq 'REML'|
   !model !eq 'reml') !THEN .
compute c2 = c .
loop l=1 to 100 .
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
*  Compute Statistics
*--------------------------------------------------------------.
compute means = T(T(T(x&*(w*p))*es)*inv(T(x&*(w*p))*x)) .
compute grpns = diag(T(x)*x) .
compute q =  T(x&*(w*p))*(es&**2)-T(T((T(x&*(w*p))*es)&**2)*inv(T(x&*(w*p))*x)) .
compute qt = csum(es&**2&*w) - csum(es&*w)**2/csum(w) .
compute qw = csum(q) .
compute qb = qt - qw .
compute dfb = ncol(x) - 1 .
compute dfw = nrow(es) - ncol(x) .
compute dft = nrow(es) - 1 .
compute se = sqrt(diag(inv(T(x&*(w*p))*x))) .
compute zvalues = means&/se .
compute pz = (1 - cdfnorm(abs(zvalues)))*2 .
compute lmeans = means - se*1.96 .
compute umeans = means + se*1.96 .
compute gmean = csum(es&*w)/csum(w) .
compute segmean = sqrt(1/csum(w)) .
compute pq = make(ncol(x),1,-9) .
loop i = 1 to ncol(x) .
+   compute pq(i,1) = 1 - chicdf(q(i,1),grpns(i,1)-1) .
end loop .

*--------------------------------------------------------------
*  Create results matrices
*--------------------------------------------------------------.
compute qtable = make(3,3,-999) .
compute qtable(1,1) = qb .
compute qtable(2,1) = qw .
compute qtable(3,1) = qt .
compute qtable(1,2) = dfb .
compute qtable(2,2) = dfw .
compute qtable(3,2) = dft .
compute qtable(1,3) = 1 - chicdf(qb,dfb) .
compute qtable(2,3) = 1 - chicdf(qw,dfw) .
compute qtable(3,3) = 1 - chicdf(qt,dft) .

compute ttable = make(1,7,-9) .
compute ttable(1,1) = gmean .
compute ttable(1,2) = segmean .
compute ttable(1,3) = gmean - segmean*1.96 .
compute ttable(1,4) = gmean + segmean*1.96 .
compute ttable(1,5) = gmean/segmean .
compute ttable(1,6) = (1 - cdfnorm(abs(gmean/segmean)))*2 .
compute ttable(1,7) = nrow(es) .

!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
compute ttable(1,1) = exp(gmean) .
compute ttable(1,3) = exp(gmean - segmean*1.96) .
compute ttable(1,4) = exp(gmean + segmean*1.96) .
!IFEND .

!IF (!print !eq 'IVZR'|!print !eq 'ivzr'|!print !eq 'Ivzr' |!print !eq 'IvZr') !THEN .
compute ttable(1,1) = (exp(2*gmean)-1)/(exp(2*gmean)+1) .
compute ttable(1,2) = -9.9999 .
compute ttable(1,3) = (exp(2*(gmean - segmean*1.96))-1)/(exp(2*(gmean - segmean*1.96))+1) .
compute ttable(1,4) = (exp(2*(gmean + segmean*1.96))-1)/(exp(2*(gmean + segmean*1.96))+1) .
!IFEND .

compute mtable = make(ncol(x),8,-9) .
compute mtable(1:ncol(x),1) = group .
compute mtable(1:ncol(x),2) = means .
compute mtable(1:ncol(x),3) = se .
compute mtable(1:ncol(x),4) = lmeans .
compute mtable(1:ncol(x),5) = umeans .
compute mtable(1:ncol(x),6) = zvalues .
compute mtable(1:ncol(x),7) = pz .
compute mtable(1:ncol(x),8) = grpns .

!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
compute mtable(1:ncol(x),2) = exp(means) .
compute mtable(1:ncol(x),4) = exp(lmeans) .
compute mtable(1:ncol(x),5) = exp(umeans) .
!IFEND .

!IF (!print !eq 'IVZR'|!print !eq 'ivzr'|!print !eq 'Ivzr' |!print !eq 'IvZr') !THEN .
compute mtable(1:ncol(x),2) = (exp(2*means)-1)/(exp(2*means)+1) .
compute mtable(1:ncol(x),3) = -9.9999*T(p) .
compute mtable(1:ncol(x),4) = (exp(2*lmeans)-1)/(exp(2*lmeans)+1) .
compute mtable(1:ncol(x),5) = (exp(2*umeans)-1)/(exp(2*umeans)+1) .
!IFEND .

compute qwtable = make(ncol(x),4,-99) .
compute qwtable(1:ncol(x),1) = group .
compute qwtable(1:ncol(x),2) = q .
compute qwtable(1:ncol(x),3) = grpns - 1 .
compute qwtable(1:ncol(x),4) = pq .

*--------------------------------------------------------------
*  Print Results
*--------------------------------------------------------------.
print /title "Version 2005.05.23".
print
  /title " *****  Inverse Variance Weighted Oneway ANOVA  *****" .
!IF (!model !eq 'ML' | !model !eq 'ml'|!model !eq 'MM'|
   !model !eq 'mm' |!model !eq 'REML' | !model !eq 'reml') !THEN .
print /title " *****  Mixed Effects Model  *****" .
!ELSE
print /title " *****  Fixed Effects Model via OLS  *****" .
!IFEND .

print qtable /format f12.4
  /title  "------- Analog ANOVA table (Homogeneity Q)  -------"
  /clabel "Q" "df" "p"
  /rlabel "Between" "Within" "Total" .
print qwtable /format f8.4
  /title  "------- Q by Group -------"
  /clabel "Group" "Qw" "df" "p" .
print ttable /format f8.4 
  /title  "------- Effect Size Results Total    -------"
  /clabel "Mean ES" "SE" "-95%CI" "+95%CI" "Z" "P" "k" 
  /rlabel "   Total" .
print mtable /format f8.4
  /title  "------- Effect Size Results by Group -------"
  /clabel "Group" "Mean ES" "SE" "-95%CI" "+95%CI" "Z" "P" "k" .

!IF (!model !eq 'MM'|!model !eq 'mm') !THEN .
 print c
   /title  '------- Method of Moments Random Effects'+
  ' Variance Component ------- '
   /rlabel "v      ="
   /format f8.5 .
!IFEND .
!IF (!model !eq 'ML' | !model !eq 'ml') !THEN .
+ compute mlc      = make(2,1,-999) .
+ compute mlc(1,1) = c .
+ compute mlc(2,1) = se_c .
+ print mlc
   /title   '------- Maximum Likelihood Random Effects'+
  ' Variance Component ------- '
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

!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
print 
  /title 'Mean ESs and 95% CIs are the exponent of the'+
  ' computed values (Odds-Ratios).' .
!IFEND .
!IF (!print !eq 'IVZR'|!print !eq 'ivzr'|!print !eq 'Ivzr' |
   !print !eq 'IvZr') !THEN .
print 
  /title 'Mean ESs and 95% CIs are the inverse Fisher'+
  ' Zr of the computed values (r).' .
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
