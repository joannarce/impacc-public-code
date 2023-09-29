%let filepath=\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\;
%include "\\rc-fs.tch.harvard.edu\crc\public\crc general information\staff folders\milliren_carly\Ozonoff - IMPACC\Reports\IMPACC Formats.sas";
**ACUTE PHASE DATA;
proc import
	datafile="&filepath\data\IMPACC_analytic_ds_2022_02_14.xlsx"
	out=acute
	dbms=xlsx replace;
run;
**CONVALESCENT SX;
proc import
	datafile="&filepath\Convalescent Paper\IMPACC_symptoms_05192022.xlsx"
	out=conv_sx
	dbms=xlsx replace;
run;
data conv_sx2;
	set conv_sx;
	studyid=pid;
	if VisitTypeID in (10,20) then delete;
	if symptom_yes=99 then symptom_yes=.;
	if symptom_yes=. then symptom_yes=0;
	**Recode no sx to 0's;
	array sx (*) loss_of_taste--organ_sum;
		do i=1 to dim(sx);
			if symptom_yes=0 & sx(i)=. then sx(i)=0;
		end;
	anosmia=max(of loss_of_taste,loss_of_smell);
	drop i;

	conv_sx_neuro=max(of anosmia,headache);
/*	conv_sx_resp=max(of red_eyes,sob,cough,sore_throat);*/
	conv_sx_cardio=max(of sob,cough);
	conv_sx_uresp=max(of red_eyes,sore_throat);
	conv_sx_sys=max(of fatigue,body_aches,fever,chills);
	conv_sx_gi=vomiting;

	sx_count=sum(of anosmia,headache,red_eyes,sob,cough,sore_throat,fatigue,body_aches,fever,chills,vomiting);
	system_count=sum(of conv_sx_neuro,conv_sx_cardio,conv_sx_uresp,conv_sx_sys,conv_sx_gi);
run;

**freq reporting organ class;
proc freq noprint
	data=conv_sx2;
	where conv_sx_neuro=1;
	table studyid /out=n_neuro;
run;
proc freq noprint
	data=conv_sx2;
	where conv_sx_uresp=1;
	table studyid /out=n_resp;
run;
proc freq noprint
	data=conv_sx2;
	where conv_sx_cardio=1;
	table studyid /out=n_cardio;
proc freq noprint
	data=conv_sx2;
	where conv_sx_sys=1;
	table studyid /out=n_sys;
run;
proc freq noprint
	data=conv_sx2;
	where conv_sx_gi=1;
	table studyid /out=n_gi;
run;
data n_system;
	merge 	n_neuro	(drop=percent rename=(count=n_neuro)) 
			n_sys 	(drop=percent rename=(count=n_sys)) 
			n_resp 	(drop=percent rename=(count=n_resp)) 
			n_cardio (drop=percent rename=(count=n_cardio))
			n_gi	(drop=percent rename=(count=n_gi)) 
	;
	by studyid;
	array n_ (*) n_neuro n_sys n_resp n_gi n_cardio;
	array g1 (*) g1_neuro g1_sys g1_resp g1_gi g1_cardio;
		do i=1 to dim(n_);
			if n_(i)>1 then g1(i)=1;	
				else if n_(i)=1 then g1(i)=0;
		end;
	drop i;
run;
**average body systems & #sx;
proc summary
	data=conv_sx2;
	class studyid;
	var sx_count system_count;
	output out=sx_avg mean=mean_sx_count mean_system_count median=median_sx_count median_system_count max=max_sx_count max_system_count;
run;
**total symptoms reported through convalesence;
proc summary
	data=conv_sx2;
	class studyid;
	var sx_count;
	output out=sx_count sum=;
run;
**ever individual sx or system;
proc summary
	data=conv_sx2;
	class studyid;
	var conv_sx_neuro
		conv_sx_gi
		conv_sx_sys
		conv_sx_uresp
		conv_sx_cardio
		symptom_yes
		anosmia	
		red_eyes		
		sob				
		cough			
		headache		
		body_aches		
		sore_throat		
		fatigue			
		fever			
		vomiting		
		chills			
	;
	output out=sx_summ max=;
run;
data all_sx;
	merge sx_summ sx_count sx_avg n_system;
	by studyid;
	if _type_=0 then delete;
	drop _type_ _freq_;

	if symptom_yes=. then delete;

	rename symptom_yes=conv_sx_any;

run;
data acute_sx;
	set acute;
	symptom_anycough=max(of symptom_cough,symptom_cough_prod);
	acute_sx_neuro=max(of symptom_anosmia,symptom_headache);
 	acute_sx_uresp=max(of symptom_conjunctivitis,symptom_sorethroat);
	acute_sx_cardio=max(of symptom_dyspnea,symptom_anycough);
	acute_sx_sys=max(of symptom_fatigue,symptom_myalgia,symptom_fever,symptom_chills);
	acute_sx_gi=symptom_nausea;
run;
**Merge all;
data all;
	merge acute_sx (in=a) all_sx (in=b);
	by studyid;
	if a;
	if b then conv=1;
		else conv=0;

	if conv=1;

	**Persistent or new sx;
	array acute(*) acute_sx_uresp acute_sx_cardio acute_sx_sys acute_sx_neuro acute_sx_gi;
	array conv_(*) conv_sx_uresp conv_sx_cardio conv_sx_sys conv_sx_neuro conv_sx_gi;
	array new (*) new_uresp new_cardio new_sys new_neuro new_gi;
	array pers(*) pers_uresp pers_cardio pers_sys pers_neuro pers_gi;
		do i=1 to dim(acute);
			if acute(i)=0 & conv_(i)=1 then new(i)=1;
				else new(i)=0;
			if acute(i)=1 & conv_(i)=1 then pers(i)=1;
				else pers(i)=0;
		end;

	**Age>=65;
	if enrollage>=65 then enrollage_c=1;
		else if .<enrollage<65 then enrollage_c=0;
		else if enrollage=. then enrollage_c=999;
	format enrollage_c enrollage_cf.;

	format sex genderf.  race racef. ethnicity ethnicityf.;

	**Combine race and ethnicity;
	raceeth=race;
	if ethnicity=1 then raceeth=6;
	format raceeth raceethf.;

	**Collapse r/e;
	raceeth2=raceeth;
	if raceeth in(1,4,9,98) then raceeth2=98;
		else if raceeth in(99,999) then raceeth2=998;
	format raceeth2 raceeth2f.;

	**Collapse race into White, Black, Other;
	if race=5 then race2=1;
		else if race=3 then race2=2;
		else if race in(1,2,4,9,98) then race2=3;
		else if race in(99,999) then race2=4;
	format race2 race2f.;

	**Collapse original race/ethnicity unknown categories;
	if race in(99,999) then race=99;
	if ethnicity in(99,999) then ethnicity=99;

	**Collapse underweight w/ normal;
	if bmicat=1 then bmicat=2;

	format baseline_img_infil imginf. bmicat bmicatf.;
run;
/*data naresh;*/
/*	set all;*/
/*	keep studyid conv_sx_any;*/
/*run;*/
/*proc freq data=naresh; table conv_sx_any; run;*/
/*proc export data=naresh outfile="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\Convalescent Stuff for acute paper\Conv Sx Data for Naresh.csv" dbms=csv replace; run;*/
**Demographic diffs in who reports during conv;
proc freq data=all; table conv_sx_any; run;
%macro freqs_any;
proc freq data=all; where &where; table &var /nocum; run;
proc freq data=all; where &where; table &var*conv_sx_any/norow nopercent chisq; run;
%mend;
%let where=							;
%let var=enrollage_c				;	%freqs_any;
%let var=sex						;	%freqs_any;
%let var=ethnicity					;	%freqs_any;
%let where=ethnicity ne 99			;	
%let var=ethnicity					;	%freqs_any;
%let var=race2						;	%freqs_any;
%let where=race2 ne 4;
%let var=race2						;	%freqs_any;
%let where=							;
%let var=comorb_htn					;	%freqs_any;
%let var=comorb_isaric_dm			;	%freqs_any;
%let var=comorb_anyresp_NOasthma	;	%freqs_any;
%let var=comorb_isaric_asthma		;	%freqs_any;	
%let var=comorb_isaric_cardiac 		;	%freqs_any;	
%let var=comorb_isaric_ckd 			;	%freqs_any;	
%let var=comorb_eversmkvape			;	%freqs_any;
%let var=comorb_countcat			;	%freqs_any;
%let var=bmicat						;	%freqs_any;
%let where=bmicat ne 999			;
%let var=bmicat						;	%freqs_any;
%let where=baseline_img_infil ne 999;
%let var=baseline_img_infil			;	%freqs_any;
%let where=							;
%let var=cluster5					;	%freqs_any;
%let var=baseline_crp_abn			;	%freqs_any;
%let var=baseline_ddimer_abn		;	%freqs_any;

**Overall symptom info;
proc freq 
	data=all;
	table conv_sx_any /nocum;
run;
title1 'Total symptoms (includes none)';
proc means n median q1 q3
	data=all;
	var median_sx_count median_system_count;
run;
title1 'Total symptoms (among those reporting symptoms)';
proc means n median q1 q3 min max
	data=all;
	where conv_sx_any=1;
	var median_sx_count median_system_count;
run;
title1 'Total symptoms (among those reporting symptoms)';
proc means n median q1 q3 min max
	data=all;
	where conv_sx_any=1;
	var max_sx_count max_system_count;
run;
title1 'Upper Respiratory';
proc freq
	data=all;
	where conv_sx_any=1;
	table (conv_sx_uresp sore_throat red_eyes)/nocum;
run;
title1 'Cardiopulmonary';
proc freq
	data=all;
	where conv_sx_any=1;
	table (conv_sx_cardio cough sob)/nocum;
run;
title1 'Systemic';
proc freq
	data=all;
	where conv_sx_any=1;
	table (conv_sx_sys fever chills fatigue body_aches)/nocum;
run;
title1 'Neuro';
proc freq
	data=all;
	where conv_sx_any=1;
	table (conv_sx_neuro headache anosmia)/nocum;
run;
title1 'GI (just vomiting)';
proc freq
	data=all;
	where conv_sx_any=1;
	table (conv_sx_gi)/nocum;
run;
title1 'Count by body system';
proc means n median q1 q3
	data=all;
	var n_resp n_cardio n_sys n_neuro n_gi;
run;
proc freq 
	data=all;
	table g1_resp g1_cardio g1_sys g1_neuro g1_gi /nocum;
run;
proc freq 
	data=all;
	table n_resp n_cardio n_sys n_neuro n_gi /nocum;
run;



**Convalescent symptoms by presenting symptoms;
%macro freq;
	proc freq data=all;
		table &acute*&conv /norow nocol agree;
	run;
%mend;
%let acute=acute_sx_uresp			;	%let conv=conv_sx_uresp		;	%freq;
%let acute=symptom_sorethroat		;	%let conv=sore_throat		;	%freq;
%let acute=symptom_conjunctivitis	;	%let conv=red_eyes			;	%freq;
%let acute=acute_sx_cardio			;	%let conv=conv_sx_cardio	;	%freq;
%let acute=symptom_anycough			;	%let conv=cough				;	%freq;
%let acute=symptom_dyspnea			;	%let conv=sob				;	%freq;
%let acute=acute_sx_sys				;	%let conv=conv_sx_sys		;	%freq;
%let acute=symptom_fever			;	%let conv=fever				;	%freq;
%let acute=symptom_chills			;	%let conv=chills			;	%freq;
%let acute=symptom_fatigue			;	%let conv=fatigue			;	%freq;
%let acute=symptom_myalgia			;	%let conv=body_aches		;	%freq;
%let acute=acute_sx_neuro			;	%let conv=conv_sx_neuro		;	%freq;
%let acute=symptom_headache			;	%let conv=headache			;	%freq;
%let acute=symptom_anosmia			;	%let conv=anosmia			;	%freq;
%let acute=acute_sx_gi				;	%let conv=conv_sx_gi		;	%freq;

proc freq data=all;
	table 	new_uresp new_cardio new_sys new_neuro new_gi
			pers_uresp pers_cardio pers_sys pers_neuro pers_gi /nocum;
run;
%macro sx_report;
/*ods select ChiSq;*/
proc freq data=all; where &where; table &var*&sx/nopercent chisq; run;
%mend;
%let sx=new_uresp					; 
%let sx=pers_uresp					;
%let sx=new_cardio					;
%let sx=pers_cardio					;
/*%let sx=new_sys 					;*/
/*%let sx=new_neuro 					;*/
/*%let sx=pers_gi						;*/
%let where=							;
%let var=enrollage_c				;	%sx_report;
%let var=sex						;	%sx_report;
%let where=ethnicity ne 99			;	
%let var=ethnicity					;	%sx_report;
%let where=race2 ne 4;
%let var=race2						;	%sx_report;
%let where=							;
%let var=comorb_htn					;	%sx_report;
%let var=comorb_isaric_dm			;	%sx_report;
%let var=comorb_anyresp_NOasthma	;	%sx_report;
%let var=comorb_isaric_asthma		;	%sx_report;	
%let var=comorb_isaric_cardiac 		;	%sx_report;	
%let var=comorb_isaric_ckd 			;	%sx_report;	
%let var=comorb_eversmkvape			;	%sx_report;
%let var=comorb_countcat			;	%sx_report;
%let where=bmicat ne 999			;
%let var=bmicat						;	%sx_report;
%let where=baseline_img_infil ne 999;
%let var=baseline_img_infil			;	%sx_report;
%let where=							;
%let var=cluster5					;	%sx_report;

**Correlation matrix of convalescent symptoms (ever);
proc corr polychoric OUTPLC=sx_corr
	data=conv_sx2;
	where symptom_yes=1;
/*	data=all;*/
/*	where conv_sx_any=1;*/
	var cough				
		sob				
		sore_throat		
		red_eyes			
		fever				
		fatigue			
		chills			
		body_aches		
		anosmia			
		headache			
		vomiting
	;
run;
proc factor
	data=sx_corr  (TYPE=CORR)
	method=principal priors=one rotate=varimax plots=scree round reorder flag=.4;
	var cough sob sore_throat red_eyes			
		fever fatigue chills body_aches 
		anosmia	headache			
		vomiting;
run;

**Inpt trajectory by convalescent symptoms;
proc freq
	data=all;
	table cluster5;
run;
title1 'Any symptoms';
proc freq 
	data=all;
	table conv_sx_any /nocum;
run;
proc freq 
	data=all;
	table conv_sx_any*cluster5/chisq nopercent norow;
run;
title1 'Total symptoms (includes none)';
proc means median q1 q3
	data=all;
	var sx_count;
run;
proc means median q1 q3
	data=all;
	class cluster5;
	var sx_count;
run;
proc npar1way kw
	data=all;
	class cluster5;
	var sx_count;
run;
title1 'Total symptoms (among those reporting symptoms)';
proc means median q1 q3
	data=all;
	where conv_sx_any=1;
	var sx_count;
run;
proc means median q1 q3
	data=all;
	where conv_sx_any=1;
	class cluster5;
	var symp_sum;
run;
proc npar1way kw
	data=all;
	where conv_sx_any=1;
	class cluster5;
	var symp_sum;
run;
title1 'Respiratory';
proc freq
	data=all;
	table (conv_sx_resp cough sore_throat sob red_eyes)/nocum;
run;
proc freq
	data=all;
	table (conv_sx_resp cough sore_throat sob red_eyes)*cluster5 /chisq nopercent norow;
run;
title1 'Systemic';
proc freq
	data=all;
	table (conv_sx_sys fever chills fatigue body_aches)/nocum;
run;
proc freq
	data=all;
	table (conv_sx_sys fever chills fatigue body_aches)*cluster5 /chisq nopercent norow;
run;
title1 'Neuro';
proc freq
	data=all;
	table (conv_sx_neuro headache anosmia)/nocum;
run;
proc freq
	data=all;
	table (conv_sx_neuro headache anosmia)*cluster5 /chisq nopercent norow;
run;
title1 'GI (just vomiting)';
proc freq
	data=all;
	table (conv_sx_gi)/nocum;
run;
proc freq
	data=all;
	table (conv_sx_gi)*cluster5 /chisq nopercent norow;
run;




**PLOT;
proc format; value visitf 0='Baseline' 3='3 months' 6='6 months' 9='9 months' 12='12 months';
value traj 1='Group 1' 2='Group 2' 3='Group 3' 4='Group 4'; run;
**Any sx by time point by trajectory group;
data conv_sx2_tr;
	merge conv_sx2 (in=a) acute (in=b keep=studyid cluster5);
	by studyid;
	if a and b;

	if visittypeid=30 then visit=3;
	if visittypeid=40 then visit=6;
	if visittypeid=50 then visit=9;
	if visittypeid=60 then visit=12;
/*	format visit visitf.;*/
	format cluster5 traj.;
run;
data acute_sx_conv_only;	
	merge acute_sx (in=a) all (in=b keep=studyid cluster5);	
	by studyid;
	if a and b;
	visit=0;
	keep  studyid visit cluster5 symptom_sorethroat	symptom_conjunctivitis	symptom_anycough			
			symptom_dyspnea	symptom_fever symptom_chills	symptom_fatigue			
			symptom_myalgia	symptom_headache symptom_anosmia symptom_nausea;
	rename 	symptom_sorethroat=sore_throat	
			symptom_conjunctivitis=red_eyes	
			symptom_anycough=cough			
			symptom_dyspnea=sob	
			symptom_fever=fever 
			symptom_chills=chills	
			symptom_fatigue=fatigue			
			symptom_myalgia=body_aches
			symptom_headache=headache 
			symptom_anosmia=anosmia 
			symptom_nausea=vomiting
	;
run;
data conv_sx2_tr2;
	set conv_sx2_tr acute_sx_conv_only;
	keep studyid visit cluster5 sore_throat red_eyes cough sob fever chills fatigue body_aches headache anosmia vomiting;
run;
proc sort data=conv_sx2_tr2; by studyid visit; run;
proc freq
	data=conv_sx2_tr;
	table cluster5*visit*Symptom_yes /out=out_sxany_traj outpct;
run;
data out_sxany_traj;
	set out_sxany_traj;
	by cluster5 visit;
	retain total;
	if first.visit then total=count;
		else total=total+count;
	y=1;
	if symptom_yes=1;
	vis_label=visit||"	(n="||trim(left(total))||")";
	pct_row_label=round(pct_row,0.1)||"%";

	rpct_row=round(pct_row/100,0.0001);
	format rpct_row percent5.1;
run;
proc template;
define statgraph sgdesign;
dynamic _VIS_LABEL _RPCT_ROW _CLUSTER5A;
dynamic _panelnumber_;
begingraph / designwidth=850 designheight=750 border=false;
   layout datalattice columnvar=_CLUSTER5A / cellwidthmin=1 cellheightmin=1 rowgutter=3 columngutter=3 rowdatarange=unionall row2datarange=unionall columndatarange=union column2datarange=unionall headerlabeldisplay=value headerlabelattrs=(size=12 family='Arial' color=CX000000) headerbackgroundcolor=white 
		columnaxisopts=( label=('Visit Month') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) discreteopts=( tickvaluefitpolicy=splitrotate)) 
		rowaxisopts=( offsetmin=0.0 label=('Any Symptom Reported (%)') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ));
      layout prototype / ;
         barchart category=_VIS_LABEL response=_RPCT_ROW / name='bar' stat=mean barwidth=0.85 barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=normal)  groupdisplay=Cluster clusterwidth=1 fillattrs=(color=CXE8E6E8 );
		endlayout;
	endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.OUT_SXANY_TRAJ template=sgdesign;
dynamic _VIS_LABEL="'VIS_LABEL'n" _RPCT_ROW="'RPCT_ROW'n" _CLUSTER5A="CLUSTER5"; 
run;

**SOB over time by traj group (Joanna);
proc freq
	data=conv_sx2_tr2;
	table cluster5*visit*sob /out=out_sxsob_traj outpct;
run;
data out_sxsob_traj;
	set out_sxsob_traj;
	by cluster5 visit;
	retain total;
	if first.visit then total=count;
		else total=total+count;
	y=1;
	if sob=1;
	vis_label=visit||"	(n="||trim(left(total))||")";
	pct_row_label=round(pct_row,0.1)||"%";

	rpct_row=round(pct_row/100,0.0001);
	format rpct_row percent5.1;
run;
proc template;
define statgraph sgdesign;
dynamic _VIS_LABEL _RPCT_ROW _CLUSTER5A;
dynamic _panelnumber_;
begingraph / designwidth=850 designheight=750 border=false;
   layout datalattice columnvar=_CLUSTER5A / cellwidthmin=1 cellheightmin=1 rowgutter=3 columngutter=3 rowdatarange=unionall row2datarange=unionall columndatarange=union column2datarange=unionall headerlabeldisplay=value headerlabelattrs=(size=12 family='Arial' color=CX000000) headerbackgroundcolor=white 
		columnaxisopts=( label=('Visit Month') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) discreteopts=( tickvaluefitpolicy=splitrotate)) 
		rowaxisopts=( offsetmin=0.0 label=('Symptom Reported (%)') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ));
      layout prototype / ;
         barchart category=_VIS_LABEL response=_RPCT_ROW / name='bar' stat=mean barwidth=0.85 barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=normal)  groupdisplay=Cluster clusterwidth=1 fillattrs=(color=CXE8E6E8 );
		endlayout;
	endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.out_sxsob_traj template=sgdesign;
where visit ne 0;
dynamic _VIS_LABEL="'VIS_LABEL'n" _RPCT_ROW="'RPCT_ROW'n" _CLUSTER5A="CLUSTER5"; 
run;



proc freq data=conv_sx2_tr2; table visit; run;
**each sx by time point (overall);
%macro sx_t;
proc freq noprint
	data=conv_sx2_tr2;
	table visit*&sx/out=sx_out outpct sparse;
run;
data sx_out;
	set sx_out;
	if &sx=1;
	length var $12;
	var="&sx";
	sys="&sys";
	drop &sx;
	keep visit count pct_row var sys;
run;
data out_sx; set out_sx sx_out; run;
%mend;
data out_sx; run;
%let sys=cardiopulm;
%let sx=cough		;	%sx_t;
%let sx=sob			;	%sx_t;	
%let sys=upperresp	;	
%let sx=sore_throat	;	%sx_t;		
%let sx=red_eyes	;	%sx_t;	
%let sys=systemic	;	
%let sx=fever		;	%sx_t;			
%let sx=fatigue		;	%sx_t;		
%let sx=chills		;	%sx_t;		
%let sx=body_aches	;	%sx_t;
%let sys=neuro		;	
%let sx=anosmia		;	%sx_t;		
%let sx=headache	;	%sx_t;	
%let sys=gi			;	
%let sx=vomiting	;	%sx_t;
proc format;
	value $varf
'cough'='Cough'		
'sob'='Shortness of breath (dyspnea)'			
'sore_throat'='Sore throat'	
'red_eyes'='Conjunctivitis/red eyes'	
'fever'='Fever'		
'fatigue'='Fatigue'		
'chills'='Chills'		
'body_aches'='Body aches (myalgia)'	
'anosmia'='Loss of smell/taste (anosmia)'		
'headache'='Headache'	
'vomiting'='Nausea/vomiting'
;	
value $sysf
'cardiopulm'='Cardiopulmonary'
'upperresp'='Upper respiratory'
'systemic'='Systemic'
'neuro'='Neurologic'
'gi'='Gastrointestinal'
;
value visit
	0='Baseline (n=589)'
	3='3 months (n=385)' 
	6='6 months (n=424)' 
	9='9 months (n=413)' 
	12='12 months (n=380)' 
;
run;
data out_sx2; 
	set out_sx;

	if var='' then delete;
	format var $varf.;

	format sys $sysf.;

	rpct_row=round(pct_row/100,0.0001);
	format rpct_row percent5.1;

	format visit visit.;

run;
proc freq data=conv_sx2_tr; table visit; run;
proc template;
define statgraph sgdesign;
dynamic _VISIT _RPCT_ROW _VAR;
begingraph / designwidth=1000 designheight=700 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none xaxisopts=( label=('Visit Month') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( viewmin=0.0 viewmax=12.0 minorticks=OFF tickvaluesequence=( start=0.0 end=12.0 increment=3.0))) 
			yaxisopts=( offsetmin=0.0 label=('Symptom Reported (%)') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ));
         seriesplot x=_VISIT y=_RPCT_ROW / group=_VAR name='series' display=(markers) connectorder=xaxis lineattrs=(pattern=SOLID thickness=2 ) markerattrs=(symbol=CIRCLEFILLED size=9 );
         discretelegend 'series' / opaque=false border=false halign=right valign=center displayclipped=true across=1 order=rowmajor location=outside valueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL);
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.OUT_SX2 template=sgdesign;
dynamic _VISIT="VISIT" _RPCT_ROW="'RPCT_ROW'n" _VAR="VAR";
run;
title;
options nodate nonumber orientation=landscape;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Convalescent Paper\Convalescent Stuff for acute paper\IMPACC Figure 3 Supplement (updated 5_19).pdf" dpi=300;
**BAR CHART;
proc template;
define statgraph sgdesign;
dynamic _VISIT _RPCT_ROW _VAR _SYS;
dynamic _panelnumber_;
begingraph / designwidth=1280 designheight=716 border=false;
   layout datapanel classvars=(_VAR) / cellwidthmin=1 cellheightmin=1 rowgutter=3 columngutter=3 rowdatarange=unionall row2datarange=unionall columndatarange=unionall column2datarange=unionall headerlabeldisplay=value rows=3 columns=4 headerlabelattrs=(size=10 family='Arial' color=CX000000) headerbackgroundcolor=white  
		columnaxisopts=( label=('Visit') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(family='Arial' size=12 style=NORMAL weight=NORMAL ) discreteopts=( tickvaluefitpolicy=splitrotate)) rowaxisopts=( offsetmin=0.0 label=('Symptom Reported (%)') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ));
      layout prototype / ;
         barchart category=_VISIT response=_RPCT_ROW / group=_SYS name='bar' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=normal) barwidth=0.85 groupdisplay=Cluster clusterwidth=0.85 /*fillattrs=(color=CXE8E6E8 )*/;
      endlayout;
    sidebar / align=bottom spacefill=false;
         discretelegend 'bar' / opaque=true border=false halign=center valign=bottom title='Organ System' displayclipped=true order=rowmajor titleattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) valueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL );
      endsidebar;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.OUT_SX2 template=sgdesign;
dynamic _VISIT="VISIT" _RPCT_ROW="'RPCT_ROW'n" _VAR="VAR" _SYS="SYS";
run;
ods pdf close;

data out_s3;
	set out_sx2;

	if sys='cardiopulm' then cardio_pct=rpct_row;
	if sys='systemic' then systemic_pct=rpct_row;
	if sys='upperresp' then upperresp_pct=rpct_row;
	if sys='gi' then gi_pct=rpct_row;
	if sys='neuro' then neuro_pct=rpct_row;
	format cardio_pct systemic_pct upperresp_pct gi_pct neuro_pct percent5.1;

	if var='cough' then cough_pct=rpct_row;	
	if var='sob' then sob_pct=rpct_row;				
	if var='sore_throat' then throat_pct=rpct_row;		
	if var='red_eyes' then eyes_pct=rpct_row;		
	if var='fever' then fever_pct=rpct_row;			
	if var='fatigue' then fatigue_pct=rpct_row;			
	if var='chills' then chills_pct=rpct_row;			
	if var='body_aches' then ache_pct=rpct_row;	
	if var='anosmia' then anosmia_pct=rpct_row;			
	if var='headache' then headache_pct=rpct_row;		
	if var='vomiting' then vomiting_pct=rpct_row;		
	format cough_pct--vomiting_pct percent5.1;
run;

proc template;
define statgraph sgdesign;
dynamic _VISIT _COUGH_PCT _VISIT2 _SOB_PCT _VISIT3 _EYES_PCT _VISIT4 _THROAT_PCT _VISIT5 _FEVER_PCT _VISIT6 _FATIGUE_PCT _VISIT7 _CHILLS_PCT _VISIT8 _ACHE_PCT _VISIT9 _HEADACHE_PCT _VISIT10 _ANOSMIA_PCT _VISIT11 _VOMITING_PCT;
begingraph / designwidth=1200 designheight=1000 border=false;
   layout lattice / rowdatarange=union columndatarange=union rows=5 columns=4 rowgutter=10 columngutter=10 rowweights=(1.0 1.0 1.0 1.0 1.0) columnweights=(1.0 1.0 1.0 1.0);
      rowaxes;
         rowaxis / display=(TICKS TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         rowaxis / display=(TICKS TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         rowaxis / display=(TICKS TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         rowaxis / display=(TICKS TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         rowaxis / display=(TICKS TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
      endrowaxes;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT response=_COUGH_PCT / name='bar' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Cough' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 style=normal weight=normal );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT2 response=_SOB_PCT / name='bar2' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Shortness of breath (dyspnea)' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay;
         entry _id='dropsite10' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay;
         entry _id='dropsite15' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay / xaxisopts=( griddisplay=off discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT3 response=_EYES_PCT / name='bar3' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Conjunctivitis/red eyes' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( griddisplay=off discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT4 response=_THROAT_PCT / name='bar4' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Sore throat' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay;
         entry _id='dropsite11' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay;
         entry _id='dropsite16' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT5 response=_FEVER_PCT / name='bar5' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Fever' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT6 response=_FATIGUE_PCT / name='bar6' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Fatigue' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT7 response=_CHILLS_PCT / name='bar7' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Chills' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT8 response=_ACHE_PCT / name='bar8' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Body aches (myalgia)' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT9 response=_HEADACHE_PCT / name='bar9' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Headache' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT10 response=_ANOSMIA_PCT / name='bar10' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Loss of  smell/taste (anosmia)' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay;
         entry _id='dropsite13' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay;
         entry _id='dropsite18' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate));
         barchart category=_VISIT11 response=_VOMITING_PCT / name='bar11' stat=mean barlabel=true barlabelattrs=(color=CX000000 family='Arial' size=8 style=NORMAL weight=NORMAL ) barwidth=0.85 groupdisplay=Cluster clusterwidth=1.0 grouporder=data fillattrs=(color=CXC6C3C6 ) outlineattrs=(color=CX000000 pattern=SOLID thickness=1 );
         entry halign=center 'Nausea/vomiting' / valign=top location=outside textattrs=(color=CX000000 family='Arial' size=12 );
      endlayout;
      layout overlay;
         entry _id='dropsite9' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay;
         entry _id='dropsite14' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      layout overlay;
         entry _id='dropsite19' halign=center '(drop a plot here...)' / valign=center;
      endlayout;
      columnaxes;
         columnaxis / display=(TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         columnaxis / display=(TICKVALUES LABEL ) label=('                                               Visit') labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         columnaxis / display=(TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
         columnaxis / display=(TICKVALUES ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL );
      endcolumnaxes;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.OUT_S3 template=sgdesign;
dynamic _VISIT="VISIT" _COUGH_PCT="'COUGH_PCT'n" _VISIT2="VISIT" _SOB_PCT="'SOB_PCT'n" _VISIT3="VISIT" _EYES_PCT="'EYES_PCT'n" _VISIT4="VISIT" _THROAT_PCT="'THROAT_PCT'n" _VISIT5="VISIT" _FEVER_PCT="'FEVER_PCT'n" _VISIT6="VISIT" _FATIGUE_PCT="'FATIGUE_PCT'n" _VISIT7="VISIT" _CHILLS_PCT="'CHILLS_PCT'n" _VISIT8="VISIT" _ACHE_PCT="'ACHE_PCT'n" _VISIT9="VISIT" _HEADACHE_PCT="'HEADACHE_PCT'n" _VISIT10="VISIT" _ANOSMIA_PCT="'ANOSMIA_PCT'n" _VISIT11="VISIT" _VOMITING_PCT="'VOMITING_PCT'n";
run;
