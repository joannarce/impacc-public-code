**use baseline_char dataset;
proc format;
	value enrollagecatf
		1='<40 years'
		2='40-49 years'
		3='50-59 years'
		4='60-69 years'
		5='70-79 years'
		6='>=80 years'
	;
	value platelets_catf
		0='<150,000/microliter'
		1='>=150,000/microliter'
		999='Missing'
	;
	value cr_catf
		0='<2 mg/dL'
		1='>=2 mg/dL'
		999='Missing'
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
	value comorb_countcat2f
		0='None'
		1='1 or more'
	;
	value discharge28_compf
		0='Not discharged'
		1='Discharged'
		2='Censored'
		3='Died'
	;
run;
data ttd;
	set allschedvisits_withbl;

	if visitnum=1;

	**KEEP ONLY ENROLLEMNTS ON OR BEFORE 12/31/2020;
	if .<baseline_date<='31DEC2020'd;

	if timetodischarge<0 then timetodischarge=.;

	**create censored indicator for discharge;
		**Not discharged;
	if timetodischarge=. then do;
		discharge28=0;
		timedischarge28=28;
	end;
		**Discharged within 28 days;
	if .<timetodischarge<=28 then do;
		discharge28=1;
		timedischarge28=timetodischarge;
	end;
		**Discharged after 28 days;
	if timetodischarge>28 then do;
		discharge28=0;
		timedischarge28=28;
	end;

	**create censored indicator for deaths;
	if diedever=. then diedever=0;
		**died within 28 days;
	if diedever=1 & timetodeath<=28 then do;
		died28=1;
		timedeath28=timetodeath;
	end;
		**died after 28 days;
	if diedever=1 & timetodeath>28 then do;
		died28=0;
		timedeath28=28;
	end;
	if diedever=0 then do;
		died28=0;
		timedeath28=28;
	end;

		**withdrew prior to 28 days;
	if .<timetowd<28 then do;
		**Censor death;
		died28=2;
		timedeath28=timetowd;
		**Censor discharge - if withdrew & not discharged within 28d;
		if discharge28=0 then do;
			discharge28=2;
			timedischarge28=timetowd;
		end;	
	end;

	**Combined competing event outcome;
	discharge28_comp=discharge28;
	timedischarge28_comp=timedischarge28;
		**death competing event;
	if discharge28=0 & died28=1 then do;
		discharge28_comp=3;
		timedischarge28_comp=timedeath28;
	end;
		**discharged then died - count as died;
	if discharge28=1 & died28=1 then do;
		discharge28_comp=3;
		timedischarge28_comp=timedeath28;
	end;
	format discharge28_comp discharge28_compf.;

	**log transform labs;
	array labs(*) baseline_lab_platelets baseline_lab_wbc baseline_lab_cr baseline_lab_lymph;
	array llabs(*) lbaseline_lab_platelets lbaseline_lab_wbc lbaseline_lab_cr lbaseline_lab_lymph;
		do i=1 to dim(labs);
			llabs(i)=log(labs(i));
		end;
	drop i;

	**binary variables for platelets and creatinine;
	if .<baseline_lab_platelets<150 then baseline_platelets_cat=1;
		else if baseline_lab_platelets>=150 then baseline_platelets_cat=0;
		else if baseline_lab_platelets=. then baseline_platelets_cat=999;
	format baseline_platelets_cat platelets_catf.;

	if .<baseline_lab_cr<2 then baseline_cr_cat=0;
		else if baseline_lab_cr>=2 then baseline_cr_cat=1;
		else if baseline_lab_cr=. then baseline_cr_cat=999;
	format baseline_cr_cat cr_catf.;

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

	baseline_sofa_minus_resp=baseline_sofa_score_new-baseline_sofa_resp;

	label 	baseline_highestresp='Baseline Respiratory Status' 
			hosp_highestresp='Respiratory Status at Hospitalization'
			hosp_spfio2ratiocat2='SpFiO2 at Hospitalization'
			baseline_spfio2ratio_lowcat2='SpFiO2 at Baseline'

		 	timedischarge28_comp='Time to Discharge (28 days)'
		  	timedeath28='Time to Death (28 days)'
	;
run;
proc sort data=ttd; by studyid; run;
**merge with file that has traj groups;
data ttd_cluster;
	merge ttd tomerge;
	by studyid;
run;

**Plots;
title1 'Baseline respiratory status';
data Risk;
   baseline_highestresp=1; output;
   baseline_highestresp=2; output;
   baseline_highestresp=3; output;
   baseline_highestresp=4; output;
   baseline_highestresp=999; output;
   format baseline_highestresp HIGHESTRESPF.;
run;
title2 'Cumulative Incidence - Discharge';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class baseline_highestresp (ref='None') dag;
   where baseline_highestresp ne 999;
   model timedischarge28_comp*discharge28_comp(0,2)=baseline_highestresp / eventcode=1;
   baseline covariates=Risk out=_null_ cif=_all_;
   id dag;
run;
title2 'Survial - 28-day mortality';
proc phreg data=ttd plots(overlay=stratum)=survival covs(aggregate);
   class baseline_highestresp (ref='None') dag;
   where baseline_highestresp ne 999;
   model timedeath28*died28(0,2)=baseline_highestresp;
   id dag;
   baseline covariates=risk;
run;

************************************************;
title1 'Baseline SpFiO2';
data Risk;
   baseline_spfio2ratio_lowcat2=1; output;
   baseline_spfio2ratio_lowcat2=2; output;
   baseline_spfio2ratio_lowcat2=3; output;
   baseline_spfio2ratio_lowcat2=999; output;
   format baseline_spfio2ratio_lowcat2 SPFIO2CAT2F.;
run;
title2 'Cumulative Incidence - Discharge';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class baseline_spfio2ratio_lowcat2 (ref='315 or higher') dag;
   where baseline_spfio2ratio_lowcat2 ne 999;
   model timedischarge28_comp*discharge28_comp(0,2)=baseline_spfio2ratio_lowcat2 / eventcode=1;
   baseline covariates=Risk out=_null_ cif=_all_;
   id dag;
run;
title2 'Survial - 28-day mortality';
proc phreg data=ttd plots(overlay=stratum)=survival covs(aggregate);
   class baseline_spfio2ratio_lowcat2 (ref='315 or higher') dag;
   where baseline_spfio2ratio_lowcat2 ne 999;
   model timedeath28*died28(0,2)=baseline_spfio2ratio_lowcat2;
   id dag;
   baseline covariates=risk;
run;

************************************************;
title1 'Hospitalization respiratory status';
data Risk;
   hosp_highestresp=1; output;
   hosp_highestresp=2; output;
   hosp_highestresp=3; output;
   hosp_highestresp=4; output;
   hosp_highestresp=999; output;
   format hosp_highestresp HIGHESTRESPF.;
run;
title2 'Cumulative Incidence - Discharge';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class hosp_highestresp (ref='None') dag;
   where hosp_highestresp ne 999;
   model timedischarge28_comp*discharge28_comp(0,2)=hosp_highestresp / eventcode=1;
   baseline covariates=Risk out=_null_ cif=_all_;
   id dag;
run;
title2 'Survial - 28-day mortality';
proc phreg data=ttd plots(overlay=stratum)=survival covs(aggregate);
   class hosp_highestresp (ref='None') dag;
   where hosp_highestresp ne 999;
   model timedeath28*died28(0,2)=hosp_highestresp;
   id dag;
   baseline covariates=risk;
run;

************************************************;
title1 'Hospitalization SpFiO2';
data Risk;
   hosp_spfio2ratiocat2=1; output;
   hosp_spfio2ratiocat2=2; output;
   hosp_spfio2ratiocat2=3; output;
   hosp_spfio2ratiocat2=999; output;
   format hosp_spfio2ratiocat2 SPFIO2CAT2F.;
run;
title2 'Cumulative Incidence - Discharge';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class hosp_spfio2ratiocat2 (ref='315 or higher') dag;
   where hosp_spfio2ratiocat2 ne 999;
   model timedischarge28_comp*discharge28_comp(0,2)=hosp_spfio2ratiocat2 / eventcode=1;
   baseline covariates=Risk out=_null_ cif=_all_;
   id dag;
run;
title2 'Survial - 28-day mortality';
proc phreg data=ttd plots(overlay=stratum)=survival covs(aggregate);
   class hosp_spfio2ratiocat2 (ref='315 or higher') dag;
   where hosp_spfio2ratiocat2 ne 999;
   model timedeath28*died28(0,2)=hosp_spfio2ratiocat2;
   id dag;
   baseline covariates=risk;
run;

************************************************;
title1 'Trajectory Grouping';
data Risk;
   cluster_traj4=1; output;
   cluster_traj4=2; output;
   cluster_traj4=3; output;
   cluster_traj4=4; output;
   format cluster_traj4 traj4f.;
run;
title2 'Cumulative Incidence - Discharge';
proc phreg data=ttd_cluster plots(overlay=stratum)=cif covs(aggregate);
   class cluster_traj4 (ref=last) dag;
   where cluster_traj4 ne .;
   model timedischarge28_comp*discharge28_comp(0,2)=cluster_traj4 / eventcode=1;
   baseline covariates=Risk out=_null_ cif=_all_;
   id dag;
run;
title2 'Survial - 28-day mortality';
proc phreg data=ttd_cluster plots(overlay=stratum)=survival covs(aggregate);
   class cluster_traj4 (ref=last) dag;
   where cluster_traj4 ne .;
   model timedeath28*died28(0,2)=cluster_traj4;
   id dag;
   baseline covariates=risk;
run;


******************************************************************************************;
**MODELS;
**Adjusted for other covariates;
**Model predicting time to discharge, mortality as competing event;
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') raceeth2 (ref='White') 
			comorb_htn (ref='No') comorb_isaric_dm (ref='No') comorb_anyresp (ref='No') comorb_isaric_cardiac (ref='No') 
			comorb_isaric_ckd (ref='No') comorb_eversmkvape (ref='No') comorb_countcat (ref='None')
			obesecat2 (ref='No') 
			bmicat (ref='Normal weight')
			comorb_isaric_obesity
			timesymptomcat2 (ref='More than 2 weeks') 
/*			baseline_spfio2ratio_lowcat2 (ref='315 or higher') */
			baseline_highestresp (ref='None') 
			baseline_platelets_cat (ref='>=150,000/microliter') baseline_cr_cat (ref='<2 mg/dL')
			/param=ref
	;
	by sex;
	model 	timedischarge28_comp*discharge28_comp(0,2)=enrollage_c sex raceeth2 
			comorb_countcat 
			comorb_htn comorb_isaric_dm	
			comorb_anyresp comorb_isaric_cardiac comorb_isaric_ckd 		
			comorb_eversmkvape 
			obesecat2 
/*			bmicat*/
			timesymptomcat2 
/*			baseline_sofa_score_new */
			baseline_sofa_minus_resp	
			baseline_highestresp
/*			baseline_spfio2ratio_lowcat2*/
			baseline_platelets_cat baseline_cr_cat 
			/risklimits=wald 
			eventcode(fg)=1;	**proportional subdistribution hazards model;
/*			eventcode(cox)=1; 	**cause-specific hazards models;*/

		id dag;
run;
proc freq 
	data=ttd;
	where obesecat2 ne 999;
	table (enrollage_c raceeth2 sex)*obesecat2 /nopercent nocol chisq;
run;

**Model predicting time to death (discharge is not a competing event);
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') raceeth2 (ref='White') 
			comorb_htn (ref='No') comorb_isaric_dm (ref='No') comorb_anyresp (ref='No') comorb_isaric_cardiac (ref='No') 
			comorb_isaric_ckd (ref='No') comorb_eversmkvape (ref='No') comorb_countcat (ref='None')
/*			obesecat2 (ref='No') */
			bmicat (ref='Normal weight')
			timesymptomcat2 (ref='More than 2 weeks') 
/*			baseline_spfio2ratio_lowcat2 (ref='315 or higher') */
			baseline_highestresp (ref='None') 
			baseline_platelets_cat (ref='>=150,000/microliter') baseline_cr_cat (ref='<2 mg/dL')
			/param=ref
	;
	model 	timedeath28*died28(0,2)=enrollage_c sex raceeth2 
/*			comorb_countcat */
			comorb_htn comorb_isaric_dm	
			comorb_anyresp comorb_isaric_cardiac comorb_isaric_ckd 		
			comorb_eversmkvape 
/*			obesecat2 */
			bmicat
			timesymptomcat2 
/*			baseline_sofa_score_new */
			baseline_sofa_minus_resp	
			baseline_highestresp
/*			baseline_spfio2ratio_lowcat2*/
			baseline_platelets_cat baseline_cr_cat 
			/risklimits=wald;
	id dag;
run;
**metabolic syndrome;
proc freq data=ttd; table sex*comorb_metabolic; run;
**Model predicting time to death (discharge is not a competing event);
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') 
			comorb_metabolic (ref='No')
			/param=ref
	;
	where sex ne 999;
	model 	timedeath28*died28(0,2)=enrollage_c sex comorb_metabolic sex*comorb_metabolic
			/risklimits=wald;
	hazardratio comorb_metabolic/diff=ref;
	id dag;
run;
proc means data=ttd; var enrollage; class sex bmicat; run;
proc freq data=ttd; table sex*bmicat*died28; run;
**time to discharge;
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') 
			comorb_metabolic (ref='No')
			/param=ref
	;
	where sex ne 999;
	model 	timedischarge28_comp*discharge28_comp(0,2)=enrollage_c sex comorb_metabolic sex*comorb_metabolic
			/risklimits=wald eventcode(fg)=1;
	hazardratio comorb_metabolic /diff=ref;
	id dag;
run;
**Explore BMI*sex more;

**Model predicting time to death (discharge is not a competing event);
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') 
			bmicat (ref='Normal weight') obesecat2 (ref='No')
			/param=ref
	;
	where sex ne 999 & obesecat2 ne 999;
	model 	timedeath28*died28(0,2)=enrollage_c sex obesecat2 sex*obesecat2
			/risklimits=wald;
	hazardratio obesecat2 /diff=ref;
	id dag;
run;
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') 
			bmicat (ref='Normal weight') obesecat2 (ref='No')
			/param=ref
	;
	where sex ne 999 & obesecat2 ne 999;
	model 	timedeath28*died28(0,2)=enrollage_c sex obesecat2 
			/risklimits=wald;
	id dag;
run;
**time to discharge;
proc phreg
	data=ttd covs(aggregate);
	class 	enrollage_c (ref='<65 years') sex (ref='Female') 
			bmicat (ref='Normal weight') obesecat2 (ref='No')
			/param=ref
	;
	where sex ne 999 & obesecat2 ne 999;
	model 	timedischarge28_comp*discharge28_comp(0,2)=enrollage_c sex obesecat2 sex*obesecat2
			/risklimits=wald eventcode(fg)=1;
	hazardratio obesecat2 /diff=ref;
	id dag;
run;

proc freq 
	data=ttd;
	where sex ne 999;
	table sex*enrollage_c sex*comorb_isaric_dm sex*comorb_htn sex*obesecat2 /nocol nopercent chisq; 
run;
proc ttest data=ttd; where sex ne 999; class sex; var enrollage; run;
proc freq data=ttd; where sex ne 999; table sex*died28/chisq; run;

ods html dpi=300;
**GENERATE FOREST PLOTS;
%macro forest;
proc import
	datafile="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Time to Discharge Models\IMPACC Time to Discharge Forest Plots.xlsx"
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
%let title1=Discharge within 28 days;	%forest;
