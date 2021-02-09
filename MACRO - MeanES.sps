*--------------------------------------------------------------
*' Macro for SPSS/Win Version 6.1 or Higher
*' Written by David B. Wilson (dwilson@crim.umd.edu)
*' Meta-Analyzes Any Type of Effect Size 
*' To use, initialize macro with the include statement:
*' INCLUDE "[drive][path]MEANES.SPS" .
*' Syntax for macro:
*' MEANES ES=varname /W=varname /PRINT=option .
*' E.g., MEANES ES = D /W = IVWEIGHT .
*' In this example, D is the name of the effect size variable
*' and IVWEIGHT is the name of the inverse variance weight
*' variable.  Replace D and INVWEIGHT with the appropriate
*' variable names for your data set.
*' /PRINT has the options "EXP" and "IVZR".  The former
*' prints the exponent of the results (odds-ratios) and
*' the latter prints the inverse Zr transform of the 
*' results.  If the /PRINT statement is ommitted, the
*' results are printed in their raw form.
*'
*' Version 2005.05.23
*'
*--------------------------------------------------------------
preserve
set printback=off
define meanes (es=!charend('/') /w=!charend('/')
  /print = !default('RAW') !charend('/'))
preserve
set printback=off mprint=off

*--------------------------------------------------------------
* Enter matrix mode and get data from active file
*--------------------------------------------------------------
matrix 
get x /file * /variables = !es !w /missing omit

*--------------------------------------------------------------
* Compute variables needed to calculate results
*--------------------------------------------------------------
compute k = nrow(x) .
compute es = make(k,1,-99) .
compute es(1:k,1) = x(1:k,1) .
compute w = make(k,1,-99) .
compute w(1:k,1) = x(1:k,2) .
release x .

*--------------------------------------------------------------
* Compute random effect variance component and new weight
*--------------------------------------------------------------
compute c = ((csum((es&**2)&*w)-csum(es&*w)**2/csum(w))-(k-1))
    /(csum(w)-csum(w&**2)/csum(w)) .
do if (c < 0) .
.compute c = 0 .
end if .
compute w_re = 1/(c + (1/w)) .

*--------------------------------------------------------------
* Calculate summary statistics
*--------------------------------------------------------------
compute df     = k - 1 .
compute mes    = csum(es&*w)   /csum(w) .
compute mes_re = csum(es&*w_re)/csum(w_re) .
compute sem    = sqrt(1/csum(w)) .
compute semre  = sqrt(1/csum(w_re)) .
compute les    = mes -    1.95996*sem .
compute ues    = mes +    1.95996*sem .
compute les_re = mes_re - 1.95996*semre .
compute ues_re = mes_re + 1.95996*semre .
compute q      = csum((es&**2)&*w)-csum(es&*w)**2/csum(w) .
do if (df>0) .
.  compute p   = 1- chicdf(q,df) .
end if .
compute z  = mes/sem .
compute z_re  = mes_re/semre .
compute pz    = (1-cdfnorm(abs(z)))*2 .
compute pz_re = (1-cdfnorm(abs(z_re)))*2 .
compute sd = sqrt(q*csum(w)**-1) .

*--------------------------------------------------------------
* Transform Output if Requested
*--------------------------------------------------------------
!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
compute mes = exp(mes) .
compute les = exp(les) .
compute ues = exp(ues) .
compute mes_re = exp(mes_re) .
compute les_re = exp(les_re) .
compute ues_re = exp(ues_re) .
compute sem = -9.9999 .
compute semre = -9.9999 .
!IFEND .

!IF (!print !eq 'IVZR'|!print !eq 'ivzr'|!print !eq 'Ivzr' 
    |!print !eq 'IvZr') !THEN .
compute mes = (exp(2*mes)-1)/(exp(2*mes)+1) .
compute les = (exp(2*les)-1)/(exp(2*les)+1) .
compute ues = (exp(2*ues)-1)/(exp(2*ues)+1) .
compute mes_re = (exp(2*mes_re)-1)/(exp(2*mes_re)+1) .
compute les_re = (exp(2*les_re)-1)/(exp(2*les_re)+1) .
compute ues_re = (exp(2*ues_re)-1)/(exp(2*ues_re)+1) .
compute sem = -9.9999 .
compute semre = -9.9999 .
!IFEND .

*--------------------------------------------------------------
* Ceate Output Matrices
*--------------------------------------------------------------
compute table1 = make(1,4,-99) .
compute table1(1,1) = k .
compute table1(1,2) = mmin(es) .
compute table1(1,3) = mmax(es) .
compute table1(1,4) = sd .

compute table2 = make(2,6,-99) .
compute table2(1,1) = mes .
compute table2(1,2) = les .
compute table2(1,3) = ues .
compute table2(1,4) = sem .
compute table2(1,5) = z .
compute table2(1,6) = pz .
compute table2(2,1) = mes_re .
compute table2(2,2) = les_re .
compute table2(2,3) = ues_re .
compute table2(2,4) = semre .
compute table2(2,5) = z_re .
compute table2(2,6) = pz_re .

compute table3 = make(1,3,-99) .
compute table3(1,1) = q .
compute table3(1,2) = df .
compute table3(1,3) = p .

*--------------------------------------------------------------
* Print summary statistics
*--------------------------------------------------------------
print /title "Version 2005.05.23".
print /title '*****  Meta-Analytic Results  *****'.
print table1
  /title '------- Distribution Description'+
  ' ---------------------------------'
  /clabel "N" "Min ES" "Max ES" "Wghtd SD"
  /format f11.3 .
print table2
  /title '------- Fixed & Random Effects Model'+
  ' -----------------------------'
  /clabel "Mean ES" "-95%CI" "+95%CI" "SE" "Z" "P"
  /rlabel "Fixed" "Random"
  /format f9.4 .
print c 
  /title '------- Random Effects Variance Component'+
  ' ------------------------'
  /rlabel 'v    =' /format f10.6 .
print table3
  /title '------- Homogeneity Analysis'+
  ' -------------------------------------'
  /clabel "Q" "df" "p"
  /format f11.4 .
print 
  /title 'Random effects v estimated via noniterative'+
  ' method of moments.' .
!IF (!print !eq 'EXP'|!print !eq 'exp'|!print !eq 'Exp') !THEN .
print 
  /title 'Mean ES and 95% CI are the exponent of the'+
  ' computed values (Odds-Ratios).' .
!IFEND .
!IF (!print !eq 'IVZR'|!print !eq 'ivzr'|!print !eq 'Ivzr' |
    !print !eq 'IvZr') !THEN .
print 
  /title 'Mean ES and 95% CI are the inverse Fisher Zr'+
  ' of the computed values (r).' .
!IFEND .
end matrix .

*--------------------------------------------------------------
* Restore settings and exit.
*--------------------------------------------------------------
restore
!enddefine
restore
