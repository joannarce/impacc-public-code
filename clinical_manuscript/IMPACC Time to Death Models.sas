**use baseline_char dataset;
proc format;
	value obesecatf	
		0='No'
		1='Obese'
		999='Missing'
	;
	value enrollagecatf
		1='<40 years'
		2='40-49 years'
		3='50-59 years'
		4='60-69 years'
		5='70-79 years'
		6='>=80 years'
	;
	/*collapse <40 with 40-49 for group 3/4 models*/
	value enrollagecat2f
		low-2='<50 years'
		3='50-59 years'
		4='60-69 years'		
		5='70-79 years'
		6='>=80 years'
	;	
	value timesymptomcat2f
		1='14 days or less'
		4='More than 2 weeks'
		99='Symptom onset date unknown'
		999='Missing'
	;	
	value age65f
		0='<65 years'
		1='>=65 years'
		999='Missing'
	;
	value comorb_countcat2f
		0='None'
		1='1 or more'
	;
run;
data ttd;
	set baseline_char;

	**KEEP ONLY ENROLLEMNTS ON OR BEFORE 12/31/2020;
	if .<baseline_date<='31DEC2020'd;

	**create censored indicator for deaths;
	if diedever=. then diedever=0;
		**died within 28 days;
	if diedever=1 & timetodeath<=28 then do;
		died28=1;
		time28=timetodeath;
	end;
		**died after 28 days;
	if diedever=1 & timetodeath>28 then do;
		died28=0;
		time28=28;
	end;
	if diedever=0 then do;
		died28=0;
		time28=28;
	end;
		**withdrew prior to 28 days;
	if .<timetowd<28 then do;
		died28=2;
		time28=timetowd;
	end;
	**create obesity collapsed variable;
	if bmicat=4 then obesecat=1;
		else if bmicat in(1,2,3) then obesecat=0;
		else if bmicat=999 then obesecat=999;
	format obesecat obesecatf.;

	**log transform labs;
	array labs(*) baseline_lab_platelets baseline_lab_wbc baseline_lab_cr baseline_lab_lymph;
	array llabs(*) lbaseline_lab_platelets lbaseline_lab_wbc lbaseline_lab_cr lbaseline_lab_lymph;
		do i=1 to dim(labs);
			llabs(i)=log(labs(i));
		end;

	**age category;
	if enrollage=. then enrollagecat=.;
		else if enrollage<40 then enrollagecat=1;
		else if 40<=enrollage<50 then enrollagecat=2;
		else if 50<=enrollage<60 then enrollagecat=3;
		else if 60<=enrollage<70 then enrollagecat=4;
		else if 70<=enrollage<80 then enrollagecat=5;
		else if enrollage>=80 then enrollagecat=6;
	format enrollagecat enrollagecatf.;

	timesymptomcat2=timesymptomcat;
	if timesymptomcat in(1,2,3) then timesymptomcat2=1;
	format timesymptomcat2 timesymptomcat2f.;

	comorb_countcat2=comorb_countcat;
	if comorb_countcat>=1 then comorb_countcat2=1;
	format comorb_countcat2 comorb_countcat2f.;

	format enrollage_c age65f.;

	baseline_sofa_minus_resp=baseline_sofa_score_new-baseline_sofa_resp;

run;
ods html dpi=300;
**overall;
proc lifetest
	data=ttd notable plots=(s);
  	time time28*died28(0,2);
run;	
**by highest resp status at baseline;
proc lifetest
	data=ttd notable plots=(s);
	where baseline_highestresp ne 999;
  	time time28*died28(0,2);
  	strata baseline_highestresp; 
run;	
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') raceeth2 (ref='White') 
			comorb_htn (ref='No') comorb_isaric_dm (ref='No') comorb_anyresp (ref='No') comorb_isaric_cardiac (ref='No') 
			comorb_isaric_ckd (ref='No') comorb_eversmkvape (ref='No') comorb_countcat (ref='None')
			obesecat (ref='No') 
			timesymptomcat2 (ref='More than 2 weeks') 
/*			baseline_spfio2ratio_lowcat2 (ref='315 or higher') */
			baseline_highestresp (ref='None') 
			/param=ref
	;
	model 	time28*died28(0,2)=enrollage_c sex raceeth2 
			comorb_countcat comorb_htn comorb_isaric_dm	
			comorb_anyresp comorb_isaric_cardiac comorb_isaric_ckd 		
			comorb_eversmkvape obesecat timesymptomcat2 
/*			baseline_sofa_score_new */
			baseline_sofa_minus_resp	
			baseline_highestresp
/*			baseline_spfio2ratio_lowcat2*/
			lbaseline_lab_platelets lbaseline_lab_cr 
			/*lbaseline_lab_wbc lbaseline_lab_lymph*/
			/risklimits=wald
	;
		id dag;
run;
**GENERATE FOREST PLOTS;
%macro forest;
proc import
	datafile="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Time to Death Models\IMPACC Time to Death Forest Plots.xlsx"
	out=plot
	dbms=xlsx replace;
run;
data plot;
	set plot;
	order=_n_;
run;
proc sort data=plot; by descending order; run;
**Get max upper limit;
proc summary data=plot; where flagunknown ne 1; var upper; output out=xmax max=xmax; run;
data xmax; set xmax; call symput('xmax',ceil(xmax)); run;
**PLOT;
proc template;
define statgraph sgdesign;
dynamic _HR _LABEL _UPPER _LOWER;
begingraph / designwidth=1000 designheight=750 border=false;
  entrytitle halign=center "&title1";
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
	layout overlay / walldisplay=none 
		xaxisopts=( display=(TICKVALUES LINE LABEL) griddisplay=on label=('Hazard Ratio (95% CI)') labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 ) 
		linearopts=( viewmin=0 viewmax=&xmax tickvaluesequence=( start=0 end=&xmax increment=2))) 
		yaxisopts=( display=(TICKVALUES LINE ) griddisplay=off tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL ) gridattrs=(color=CX848284 pattern=2 thickness=1 ) discreteopts=( tickvaluefitpolicy=none));
         scatterplot x=_HR y=_LABEL / group=_REF xerrorupper=_UPPER xerrorlower=_LOWER name='scatter' grouporder=data markerattrs=(color=CX000000 symbol=CIRCLEFILLED size=7 ) errorbarattrs=(color=CX000000 pattern=SOLID thickness=1 );
         referenceline x=1.0 / name='vref' xaxis=X curvelabelposition=max lineattrs=(color=CX7d7b7a pattern=solid thickness=1 );
      endlayout;
   endlayout;
endgraph;
end;
run;
proc sgrender data=WORK.plot template=sgdesign;
where flagunknown ne 1;
dynamic _HR="'HR'n" _LABEL="LABEL" _UPPER="UPPER" _LOWER="LOWER";
run;
%mend;
/*ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Trajectory\IMPACC Trajectory Grouping Forest Plots.pdf";*/
%let title1=28-day Mortality;	%forest;

**diff mortality by site;
proc freq data=ttd; table dag*died28 /chisq nocol nopercent;run;
proc lifetest
	data=ttd notable plots=(s);
  	time time28*died28(0,2);
  	strata dag; 
run;	
