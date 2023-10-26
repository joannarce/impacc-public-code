**GENERATE FIGURE 2A;
**FOREST PLOTS;
title1 'Predictors of Deficit Cluster vs. not';
proc import
	datafile="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\Convalescent Logistic Regression.xlsx"
	out=figs
	dbms=xlsx replace;
	sheet="upd logistic";
run;
data figs;
	set figs;

	if UNK=1 then delete;
run;
proc sort data=figs; by descending varnum descending varlevel; run;
ods html dpi=200;
proc format;
	value logdisplayf
		-5		='0'
		-3		='0.05'
		-2.3	='0.10'
		-1.61	='0.20'
		-1.39	='0.25'
		-0.69	='0.50'
		-0.29	='0.75'
		0.00	='1'
		0.41	='1.5'
		0.69	='2'
		0.92	='2.5'
		1.10	='3'
		1.39	='4'
		1.61	='5'
		1.79 	='6'
		1.95	='7'
		2.08	='8'
		2.20	='9'
		2.30	='10'
	;
run; 
**Figure 2A;
title ' ';
options nodate nonumber orientation=landscape;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\IMPACC Conv Figure 2A.pdf" dpi=300;
proc template;
define statgraph sgdesign;
dynamic _ESTIMATE _LABEL _UPPER _LOWER;
begingraph / designwidth=1000 designheight=900 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / walldisplay=none 
		xaxisopts=( label=('Odds Ratio (95% CI)') labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=NORMAL ) linearopts=(tickvalueformat=LOGDISPLAYF. tickvaluelist=(-1.39 -0.69 -0.29 0.00 0.41 0.69 1.10 1.39))) 
		yaxisopts=( display=(LINE TICKVALUES TICKS ) griddisplay=on tickvalueattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) discreteopts=( tickvaluefitpolicy=none));
         scatterplot x=_ESTIMATE y=_LABEL / xerrorupper=_UPPER xerrorlower=_LOWER name='scatter' markerattrs=(color=CX000000 symbol=CIRCLEFILLED size=9 ) errorbarattrs=(color=CX000000 pattern=SOLID thickness=1 );
         referenceline x=0.0 / name='vref' xaxis=X curvelabelposition=max lineattrs=(color=CX848284 pattern=MEDIUMDASH thickness=1 );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.FIGS template=sgdesign;
dynamic _ESTIMATE="ESTIMATE" _LABEL="LABEL" _UPPER="UPPER" _LOWER="LOWER";
run;
ods pdf close;
