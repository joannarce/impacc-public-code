**Generate suppFigure2;
**Figure for presenting and convalescent sx;
%macro freq;
proc freq noprint
	data=all_final;
	table &presenting /out=outpresent outpct;
run;
data outpresent2;
	set outpresent;
	if &presenting=1;
	drop &presenting;
run;
proc freq noprint
	data=all_final;
	table &conv /out=outconv outpct;
run;
data outconv2;
	set outconv;
	if &conv=1;
	drop &conv;
run;
data allsx;
	length vartext $25.;
	set outpresent2 (in=a) outconv2 (in=b);
	if a then t=0;
		else if b then t=1;
	format t timef.;
	vartext="&presenting";
	varnum=&varnum;
	format vartext $vartextf.;
run;
data allsx_final; set allsx_final allsx; run;
%mend;
proc format;
	value $vartextf
		"acute_sx_uresp"			="Any upper respiratory (sore throat, conjunctivitis)"			
		"symptom_conjunctivitis"	="     Conjunctivitis/red eyes"	
		"symptom_sorethroat"		="     Sore throat"		
		"acute_sx_cardio"			="Any cardiopulmonary (cough, dyspnea)"			
		"symptom_anycough"			="     Dyspnea"		
		"symptom_dyspnea"			="     Cough"			
		"acute_sx_sys"				="Any systemic (fever, fatigue, myalgia, chills)"			
		"symptom_fatigue"			="     Fatigue"			
		"symptom_myalgia"			="     Myalgia"			
		"symptom_fever"				="     Fever"			
		"symptom_chills"			="     Chills"			
		"acute_sx_neuro"			="Any neurologic (headache, anosmia)"			
		"symptom_headache"			="     Headache"		
		"symptom_anosmia"			="     Anosmia"			
		"acute_sx_gi"				="Any gastrointestinal (nausea/vomiting)"		
	;
	value timef
		0="Presenting"
		1="Convalescent"
	;
run;
data allsx_final; run;
%let presenting=acute_sx_uresp			;	%let conv=conv_sx_uresp	;	%let varnum=1;		%freq;
%let presenting=symptom_conjunctivitis	;	%let conv=red_eyes		;	%let varnum=2;		%freq; 
%let presenting=symptom_sorethroat		;	%let conv=sore_throat	;	%let varnum=3;	 	%freq; 
%let presenting=acute_sx_cardio			;	%let conv=conv_sx_cardio;	%let varnum=4;		%freq;
%let presenting=symptom_anycough		;	%let conv=cough			;	%let varnum=5;	 	%freq; 
%let presenting=symptom_dyspnea			;	%let conv=sob			;	%let varnum=6;	 	%freq; 
%let presenting=acute_sx_sys			;	%let conv=conv_sx_sys	;	%let varnum=7;		%freq;
%let presenting=symptom_fatigue			;	%let conv=fatigue		;	%let varnum=8;	 	%freq; 
%let presenting=symptom_myalgia			;	%let conv=body_aches	;	%let varnum=9;	 	%freq; 
%let presenting=symptom_fever			;	%let conv=fever			;	%let varnum=10;	 	%freq; 
%let presenting=symptom_chills			;	%let conv=chills		;	%let varnum=11;	 	%freq;
%let presenting=acute_sx_neuro			;	%let conv=conv_sx_neuro	;	%let varnum=12;		%freq;
%let presenting=symptom_headache		;	%let conv=headache		;	%let varnum=13;	 	%freq;
%let presenting=symptom_anosmia			;	%let conv=anosmia		;	%let varnum=14;	 	%freq;
%let presenting=acute_sx_gi				;	%let conv=conv_sx_gi	;	%let varnum=15;	 	%freq;	
data allsx_final2; set allsx_final; if vartext='' then delete; run;
proc print data=allsx_final2 noobs; run;
proc sort data=allsx_final2; by descending varnum descending t; run;

ods html style=journal dpi=300; 
**Figure 2S;
title ' ';
options nodate nonumber orientation=landscape;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\IMPACC Conv Figure 2S.pdf" dpi=300 style=journal;
proc template;
define statgraph sgdesign;
dynamic _VARTEXT _PERCENT _T;
begingraph / designwidth=1000 designheight=850 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none xaxisopts=( offsetmin=0.0 label=('Percent') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL )) yaxisopts=( display=(TICKS TICKVALUES LINE ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) discreteopts=( tickvaluefitpolicy=none));
         barchart category=_VARTEXT response=_PERCENT / group=_T name='bar(h)' stat=mean orient=horizontal groupdisplay=Cluster clusterwidth=0.85;
         discretelegend 'bar(h)' /autoitemsize=true opaque=false border=false halign=right valign=top displayclipped=true down=1 order=columnmajor location=inside valueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.ALLSX_FINAL2 template=sgdesign;
dynamic _VARTEXT="VARTEXT" _PERCENT="PERCENT" _T="T";
run;
ods pdf close;
