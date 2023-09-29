proc format;
	value discharge28_compf
		0='Not discharged'
		1='Discharged'
		2='Censored'
		3='Died'
	;
	value clusterf
		1='1'
		2='2'
		3='3'
		4='4'
		5='5'
	;
run;
data ttd;
	set t;	**note this dataset is created in main cleaning file in order to calculate death by qtr;
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
		**died within 28 days;
	if died28=1 then do;
		timedeath28=timetodeath;
	end;
		**died after 28 days;
	if died28=0 & timetodeath>28 then do;
		timedeath28=28;
	end;
	if died28=0 then do;
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

	label 	timedischarge28_comp='Time to Discharge (28 days)'
		  	timedeath28='Time to Death (28 days)'
	;
run;
proc sort data=ttd; by studyid; run;
**Plots;
title1 'TG';
data Risk;
   cluster5=1; output;
   cluster5=2; output;
   cluster5=3; output;
   cluster5=4; output;
   cluster5=5; output;
   format cluster5 best12.;
run;
title2 'Cumulative Incidence - Discharge - BY TG';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class cluster5 dag;
   model timedischarge28_comp*discharge28_comp(0,2)=cluster5 / eventcode=1;
   baseline covariates=Risk out=disch_traj cif=_all_;
   label cluster5 ='Trajectory Group';
   id dag;
run;
**manipulate so each traj stored in a separate var;
data disch_traj2;
	set disch_traj;

	array cif_ (*) cif_1-cif_5;
		do i=1 to dim(cif_);
			if cluster5=i then cif_(i)=cif;
		end;
run;
proc summary 
	data=disch_traj2;
	class timedischarge28_comp;
	var cif_1-cif_5;
	output out=disch_traj3 max=;
run;
data disch_traj3; set disch_traj3; where _type_=1; drop _type_ _freq_; run;
/*ods html dpi=300;*/

options nodate nonumber;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Manuscript Draft\For journal\vector based graphs\IMPACC Suppl Fig 1D.pdf" dpi=300;
title "	";
proc template;
define statgraph sgdesign;
dynamic _TIMEDISCHARGE28_COMP _CIF_1A _TIMEDISCHARGE28_COMP2 _CIF_2A _TIMEDISCHARGE28_COMP3 _CIF_3A _TIMEDISCHARGE28_COMP4 _CIF_4A _TIMEDISCHARGE28_COMP5 _CIF_5A;
begingraph / designwidth=750 designheight=550 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none xaxisopts=( label=('Days from admission') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( viewmin=0.0 viewmax=28.0 tickvaluepriority=TRUE tickvalueformat=BEST6. tickvaluelist=(0.0 4.0 7.0 14.0 21.0 28.0))) yaxisopts=( label=('Proportion Discharged') labelattrs=(color=CX000000 family='Arial' size=14 ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( viewmin=-0.01 viewmax=1.0));
         stepplot x=_TIMEDISCHARGE28_COMP y=_CIF_1A / name='step' 	legendlabel='1' connectorder=xaxis lineattrs=(color=CX639A21 pattern=SOLID thickness=2 );
         stepplot x=_TIMEDISCHARGE28_COMP2 y=_CIF_2A / name='step2' legendlabel='2' connectorder=xaxis lineattrs=(color=CX39828C pattern=SHORTDASH thickness=2 );
         stepplot x=_TIMEDISCHARGE28_COMP3 y=_CIF_3A / name='step3' legendlabel='3' connectorder=xaxis lineattrs=(color=CX6371AD pattern=MEDIUMDASH thickness=2 );
         stepplot x=_TIMEDISCHARGE28_COMP4 y=_CIF_4A / name='step4' legendlabel='4' connectorder=xaxis lineattrs=(color=CXBD7D31 pattern=DASH thickness=2 );
         stepplot x=_TIMEDISCHARGE28_COMP5 y=_CIF_5A / name='step5' legendlabel='5' connectorder=xaxis lineattrs=(color=CX9C3418 pattern=LONGDASH thickness=2 );
         discretelegend 'step' 'step2' 'step3' 'step4' 'step5' / opaque=false border=false halign=center valign=top displayclipped=true down=1 order=columnmajor location=outside 
			title='Trajectory Group' titleattrs=(family='Arial' size=12 style=NORMAL weight=BOLD ) valueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL );
      endlayout;
   endlayout;
endgraph;
end;
run;
proc sgrender data=WORK.DISCH_TRAJ3 template=sgdesign;
dynamic _TIMEDISCHARGE28_COMP="'TIMEDISCHARGE28_COMP'n" _CIF_1A="'CIF_1'n" _TIMEDISCHARGE28_COMP2="'TIMEDISCHARGE28_COMP'n" _CIF_2A="'CIF_2'n" _TIMEDISCHARGE28_COMP3="'TIMEDISCHARGE28_COMP'n" _CIF_3A="'CIF_3'n" _TIMEDISCHARGE28_COMP4="'TIMEDISCHARGE28_COMP'n" _CIF_4A="'CIF_4'n" _TIMEDISCHARGE28_COMP5="'TIMEDISCHARGE28_COMP'n" _CIF_5A="'CIF_5'n";
run;
ods pdf close;

title2 'Cumulative Incidence - Discharge - OVERALL';
proc phreg data=ttd plots(overlay=stratum)=cif covs(aggregate);
   class dag;
   model timedischarge28_comp*discharge28_comp(0,2)= / eventcode=1;
   baseline /*covariates=Risk*/ out=disch cif=_all_;
   id dag;
run;
options nodate nonumber;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Manuscript Draft\For journal\vector based graphs\IMPACC Suppl Fig 1C.pdf" dpi=300;
title "	";
proc template;
define statgraph sgdesign;
dynamic _TIMEDISCHARGE28_COMP _CIF;
begingraph / designwidth=750 designheight=550 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none 
			xaxisopts=( display=(TICKS TICKVALUES LINE LABEL ) label=('Days from admission') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( minorticks=OFF tickvaluepriority=TRUE tickvalueformat=BEST6. tickvaluelist=(0.0 4.0 7.0 14.0 21.0 28.0))) 
			yaxisopts=( label=('Proportion Discharged') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL )  linearopts=( viewmin=0.0 viewmax=1.0 minorticks=OFF tickvaluesequence=( start=0.0 end=1.0 increment=0.2)));
         stepplot x=_TIMEDISCHARGE28_COMP y=_CIF / name='step' connectorder=xaxis lineattrs=(color=CX000000 pattern=SOLID thickness=2 );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.DISCH template=sgdesign;
dynamic _TIMEDISCHARGE28_COMP="'TIMEDISCHARGE28_COMP'n" _CIF="CIF";
run;
ods pdf close;

title2 'Survial - 28-day mortality';
proc phreg data=ttd plots(overlay=stratum)=survival covs(aggregate);
   model timedeath28*died28(0,2)=;
   id dag;
   output out=surv survival=s;
run;
options nodate nonumber;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Manuscript Draft\For journal\vector based graphs\IMPACC Suppl Fig 1A.pdf" dpi=300;
title "	";
proc template;
define statgraph sgdesign;
dynamic _TIMEDEATH28A _S;
begingraph / designwidth=750 designheight=550 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none 
			xaxisopts=( offsetmin=0.0 display=(TICKS TICKVALUES LINE LABEL ) label=('Days from admission') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( minorticks=OFF tickvaluepriority=TRUE tickvalueformat=BEST6. tickvaluelist=(0.0 4.0 7.0 14.0 21.0 28.0))) 		
			yaxisopts=( offsetmin=0.0 label=('Survival Proportion') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( viewmin=0.0 viewmax=1.0));
         stepplot x=_TIMEDEATH28A y=_S / name='step' connectorder=xaxis lineattrs=(color=CX000000 pattern=SOLID thickness=2 );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.SURV template=sgdesign;
dynamic _TIMEDEATH28A="TIMEDEATH28" _S="S";
run;
ods pdf close;


**Mortality by quarter;
**Deaths by site and calendar qtr;
data t;
	set traj_all; 
	where visitday=0;
	if cluster5=5 then died28=1;
		else died28=0;
	qtryear=baseline_date;
	format qtryear yyq.;

	if race2 in(2,3) or ethnicity=1 then nonwhite=1;
		else if race2=1 and ethnicity=0 then nonwhite=0;

	**ltf within 28 days;
	if .<timetowd<=28 then ltf_28=1;
		else ltf_28=0;
run;
**Differential loss to follow-up by site, quarter or race/eth;
proc freq data=t; table dag*ltf_28 /chisq nocol nopercent; run;
proc freq data=t; table qtryear*ltf_28 /chisq nocol nopercent; run;
proc freq data=t; table nonwhite*ltf_28 /chisq nocol nopercent; run;

**Differences by site/quarter in 28-d mortality;
proc freq data=t; table dag*died28 /chisq nocol nopercent; run;
proc freq data=t; table qtryear*died28 /chisq nopercent nocol out=qtrdeath outpct; run;
data qtrdeath2;
	set qtrdeath;
	by qtryear;

	retain tot;
	if first.qtryear then tot=count;
		else tot=tot+count;
	if died28=1;

	qtryear_n=vvalue(qtryear)||" 			(enrolled n="||trim(left(tot))||")";

	count_n="n="||trim(left(count));
run;
/*ods html dpi=300;*/
options nodate nonumber;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript\Manuscript Draft\For journal\vector based graphs\IMPACC Suppl Fig 1B_noannotation.pdf" dpi=300;
title "	";
proc template;
define statgraph sgdesign;
dynamic _PCT_ROW _PCT_ROW2 _COUNT_N _QTRYEAR_N _QTRYEAR_N2;
begingraph / designwidth=750 designheight=550 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none 
					xaxisopts=( display=(LINE TICKVALUES LABEL ) label=('Enrollment Quarter') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=11 style=NORMAL weight=NORMAL ) discreteopts=( tickvaluefitpolicy=split)) 
					yaxisopts=( offsetmin=0.0 label=('Mortality (%)') labelattrs=(color=CX000000 family='Arial' size=14 style=NORMAL weight=NORMAL ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=( viewmin=0.0 viewmax=12.0 minorticks=OFF tickvaluesequence=( start=0.0 end=12.0 increment=2.0)));
         barchart category=_QTRYEAR_N response=_PCT_ROW / name='bar' stat=mean groupdisplay=Cluster clusterwidth=1.0 fillattrs=(color=CXC6C3C6 ) outlineattrs=(pattern=SOLID thickness=1 );
         scatterplot x=_QTRYEAR_N2 y=_PCT_ROW2 / datalabel=_COUNT_N name='scatter' labelstrip=true datalabelposition=TOP markerattrs=(size=0 ) datalabelattrs=(color=CX000000 family='Arial' size=11 style=NORMAL weight=NORMAL );
/*         entry halign=right 'p=0.84' / valign=top textattrs=(color=CX000000 family='Arial' size=12 style=normal weight=normal );*/
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.QTRDEATH2 template=sgdesign;
dynamic _PCT_ROW="'PCT_ROW'n" _PCT_ROW2="'PCT_ROW'n" _COUNT_N="COUNT_N" _QTRYEAR_N="'QTRYEAR_N'n" _QTRYEAR_N2="'QTRYEAR_N'n";
run;
ods pdf close;

**is site mortality diff d/t site r/eth breakdown?;
proc freq data=t; table dag*nonwhite*died28 /chisq nocol nopercent cmh ; run;
