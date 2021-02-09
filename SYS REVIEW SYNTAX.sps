* Encoding: UTF-8.

*The below syntax was used on the RAW data set to assign value labels, rename variables, and create summed score variables.

DELETE VARIABLES StartDate
EndDate
Status
IPAddress
Progress
Duration__in_seconds_
Finished
RecordedDate
ResponseId
RecipientLastName
RecipientFirstName
RecipientEmail
ExternalReference
LocationLatitude
LocationLongitude
DistributionChannel
UserLanguage.
EXECUTE.

COMPUTE STUDYID=0.
FORMATS STUDYID(F6.0).
EXECUTE.

SORT CASES BY Q2(A) Q1(A).
EXECUTE.
**Manually number each included study starting with 1.

RENAME VARIABLES (Q1 Q2 Q3 Q4 Q5_1 Q5_2 Q5_3 Q5_4 Q5_4_TEXT Q6_1 Q6_2 Q6_3 Q6_4 Q6_5 Q6_6 Q6_7
Q6_7_TEXT Q7_1 Q7_2 Q7_3 Q7_4 Q7_5 Q7_6 Q7_7 Q7_6_TEXT Q8 Q9 Q9_6_TEXT Q10_1 Q10_2 Q10_3 Q10_4 Q10_5
Q10_6 Q10_7 Q10_8 Q10_9Q10_8_TEXT Q11 Q12_1 Q12_2 Q12_2_TEXT Q13 Q14_1.0 Q14_2.0 Q15_1 Q15_2 Q15_2_TEXT 
Q16 Q17_1.0 Q17_2.0 Q18 Q19_1 Q19_2 Q19_2_TEXT Q21_1.0 Q21_2.0 Q25 = AUTHOR YEAR N AGE OLDGENDER1 
OLDGENDER2 OLDGENDER3 OLDGENDER4 GENDEROTH OLDRACE1 OLDRACE2 OLDRACE3 OLDRACE4 OLDRACE5 
OLDRACE6 OLDRACE7 RACEOTH OLDSEXUALITY1 OLDSEXUALITY2 OLDSEXUALITY3 OLDSEXUALITY4 OLDSEXUALITY5 
OLDSEXUALITY6 OLDSEXUALITY7 SEXUALITYOTH LOCATION DESIGN DESIGNOTH OLDTHEORY1 OLDTHEORY2
OLDTHEORY3 OLDTHEORY4 OLDTHEORY5 OLDTHEORY6 OLDTHEORY7 OLDTHEORY8 OLDTHEORY9 THEORYOTH SKILLASSESS 
OLDOBVINS1 OLDOBSVIN2 OBVINSTROTH OBVSTEPS OBVREL OBVVAL OLDPROXINSTR1 OLDPROXINSTR2 PROXINSTROTH 
PROXSTEPS PROXREL PROXVAL EFFOPER OLDEFFINSTR1 OLDEFFINSTR2 EFFINSTROTH EFFREL EFFVAL RESULTS).
EXECUTE.

**Compute new variables and values for all of the multiple choice answers.
COMPUTE GENDER=0.
COMPUTE GENDERSUM=SUM(Q5_1 to Q5_4).
IF (Q5_1=1) GENDER=1.
IF (Q5_2=1) GENDER=2.
IF (Q5_3=1) GENDER=3.
IF (Q5_4=1) GENDER=4.
IF (GENDERSUM GT 1) GENDER=5.
EXECUTE.

VALUE LABELS 
GENDER
1 'Male'
2 'Female'
3 'Transgender'
4 'Other'
5 'Multiple genders'.
EXECUTE.

COMPUTE RACE=0.
COMPUTE RACESUM=SUM(Q6_1 to Q6_7).
IF (Q6_1=1) RACE=1.
IF (Q6_2=1) RACE=2.
IF (Q6_3=1) RACE=3.
IF (Q6_4=1) RACE=4.
IF (Q6_5=1) RACE=5.
IF (Q6_6=1) RACE=6.
IF (Q6_7=1) RACE=7.
IF (RACESUM GT 1) RACE=8.
EXECUTE.

VALUE LABELS 
RACE
1 'American Indian/Alaska Native'
2 'Asian'
3 'Black/African-American'
4 'Hawaiian/Other Pacific Islander'
5 'White'
6 'Hispanic/Latino'
7 'Other'
8 'Multiple races'
0 'Not Stated'.
EXECUTE.

COMPUTE SEXUALITY=0.
COMPUTE SEXUALITYSUM=SUM(Q7_1 to Q7_6).
IF (Q7_1=1) SEXUALITY=1.
IF (Q7_2=1) SEXUALITY=2.
IF (Q7_3=1) SEXUALITY=3.
IF (Q7_4=1) SEXUALITY=4.
IF (Q7_5=1) SEXUALITY=5.
IF (Q7_6=1) SEXUALITY=6.
IF (Q7_7=1) SEXUALITY=7.
IF (SEXUALITYSUM GT 1) SEXUALITY=8.
EXECUTE.

VALUE LABELS 
SEXUALITY
1 'Heterosexual'
2 'Homosexual'
3 'Bisexual'
4 'Asexual'
5 'Transexual'
6 'Other'
7 'Not Stated'
8 'Multiple sexualities'.
EXECUTE.

COMPUTE THEORY=0.
COMPUTE THEORYSUM=SUM(Q10_1 to Q10_8).
IF (Q10_1=1) THEORY=1.
IF (Q10_2=1) THEORY=2.
IF (Q10_3=1) THEORY=3.
IF (Q10_4=1) THEORY=4.
IF (Q10_5=1) THEORY=5.
IF (Q10_6=1) THEORY=6.
IF (Q10_7=1) THEORY=7.
IF (Q10_8=1) THEORY=8.
IF (Q10_9=1) THEORY=9.
IF (THEORYSUM GT 1) THEORY=10.
EXECUTE.

VALUE LABELS 
THEORY
1 'Health Belief Model'
2 'Social Cognitive Theory/Self-Efficacy Theory'
3 'Transtheoretical Model'
4 'Social Ecological Model'
5 'Diffusion of Innovation'
6 'Theory of Planned Behavior'
7 'Theory of Reasoned Action'
8 'Other'
9 'Not Stated'
10 'Multiple theories'.
EXECUTE.

COMPUTE OBVINSTR=0.
IF (Q12_1=1) OBVINSTR=1.
IF (Q12_2=1) OBVINSTR=2.
EXECUTE.

VALUE LABELS
OBVINSTR
1 'Measure of Observed Condom Use Skills'
2 'Other'
0 'No direct observation'.
EXECUTE.

COMPUTE PROXINSTR=0.
IF (Q15_1=1) PROXINSTR=1.
IF (Q15_2=1) PROXINSTR=2.
EXECUTE.

VALUE LABELS
PROXINSTR
1 'Condom Use Skills Checklist (CUSC)'
2 'Other'
0 'No proxy instrument'.
EXECUTE.

COMPUTE EFFINSTR=0.
IF (Q19_1=1) EFFINSTR=1.
IF (Q19_2=1) EFFINSTR=2.
EXECUTE.

VALUE LABELS
EFFINSTR
1 'Condom Use Self-Efficacy Scale (CUSES)'
2 'Other'.
EXECUTE.

VARIABLE LABELS
STUDYID 'ID for this publication'
AUTHOR 'First Author Last Name'
YEAR 'Publication Year'
N 'Sample N'
AGE 'Age Range or Mean (SD)'
GENDER 'Gender categories included in the sample'
RACE 'Race/Ethnicity categories included in sample'
SEXUALITY 'Sexual orientation categories included in sample'
LOCATION 'Study Location'
DESIGN 'Study Design'
THEORY 'Theories or Frameworks used in the study'
SKILLASSESS 'How were skills assessed?'
OBVINSTR 'What direct observation instrument was used?'
OBVSTEPS 'Number of steps in observation instrument'
OBVREL 'Observation Instrument Reliability'
OBVVAL 'Observation Instrument Validity'
PROXINSTR 'What proxy instrument was used?'
PROXSTEPS 'Number of steps in proxy instrument'
PROXREL 'Proxy Instrument Reliability'
PROXVAL 'Proxy Instrument Validity'
EFFOPER 'How was self-efficacy operationalized?'
EFFINSTR 'What instrument was used to assess this self-efficacy?'
EFFREL 'Self-efficacy Instrument Reliability'
EFFVAL 'Self-efficacy Instrument Validity'
RESULTS 'Results reported'
CORREFF1 'First reported correlation coefficient'
CORREFF2 'Second reported correlation coefficient'
CORREFF3 'Third reported correlation coefficient'
GENDERSUM 'Total number of distinct genders in sample'
RACESUM 'Total number of distinct races/ethnicities in sample'
SEXUALITYSUM 'Total number of distinct sexualities in sample'
THEORYSUM 'Total number of distinct theories used in study'.
EXECUTE.

COMPUTE CORREFF1=0.
COMPUTE CORREFF2=0.
COMPUTE CORREFF3=0.
FORMATS CORREFF1 to CORREFF3 (F6.3).
EXECUTE.
*Manually extract the correlation coeeficients from the results into a new dataset that will contain only studies with valid effect sizes for meta-analysis.



