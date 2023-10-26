**GENERATE FIGURE 2B;
**PLOT Multinomial;
**FOREST PLOTS;
proc format; 
	value compf 2='Physical predominant (2) vs. minimal deficit (1)' 
				3='Mental/cognitive predominant (5) vs. minimal deficit (1)' 
				4='Multi/pan domain (6) vs. minimal deficit (1)'; 
	value compnewf
				2='PHY vs. MIN'
				3='COG vs. MIN'
				4='MLT vs. MIN'
	;
run;
title1 'Predictors of individual Deficit Clusters vs. minimal/none';
proc import
	datafile="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\Convalescent Logistic Regression.xlsx"
	out=multi
	dbms=xlsx replace;
	sheet="upd multi";
run;
data multi;
	set multi;

	if UNK=1 then delete;
	comp=cluster4;
	format comp compnewf.;
run;
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
proc sort data=multi; by comp descending varnum descending varlevel; run;
**Figure 2B;
title ' ';
options nodate nonumber orientation=landscape;
ods noproctitle;
ods graphics on /outputfmt=pdf; 
ods pdf file="Z:\Ozonoff - IMPACC\Manuscript\Convalescent Paper\IMPACC Conv Figure 2B.pdf" dpi=300;
proc template;
define statgraph sgdesign;
dynamic _COMP _ESTIMATE _LABEL _UPPER _LOWER;
dynamic _panelnumber_;
begingraph / designwidth=1400 designheight=900 border=false;
 layout datalattice columnvar=_COMP / cellwidthmin=1 cellheightmin=1 rowgutter=3 columngutter=3 rowdatarange=unionall row2datarange=unionall columndatarange=union column2datarange=union headerlabeldisplay=value headerlabelattrs=(size=10 family='Arial' color=CX000000) headerbackgroundcolor=white
		columnaxisopts=( label=('Odds Ratio (95% CI)') labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL ) linearopts=(tickvaluefitpolicy=none tickvalueformat=LOGDISPLAYF. tickvaluelist=(-2.3 -1.39 -0.69 0.00 0.69 1.10 1.39 1.79 2.08 2.30))) 
		rowaxisopts=( display=(LINE TICKVALUES TICKS ) griddisplay=on tickvalueattrs=(color=CX000000 family='Arial' size=9 style=NORMAL weight=BOLD ) discreteopts=( tickvaluefitpolicy=none));
      layout prototype /;          
		scatterplot x=_ESTIMATE y=_LABEL / xerrorupper=_UPPER xerrorlower=_LOWER name='scatter' markerattrs=(color=CX000000 symbol=CIRCLEFILLED size=9 ) errorbarattrs=(color=CX000000 pattern=SOLID thickness=1 );
        referenceline x=0.0 / name='vref' xaxis=X curvelabelposition=max lineattrs=(color=CX848284 pattern=MEDIUMDASH thickness=1 );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.multi template=sgdesign;
dynamic _COMP='COMP' _ESTIMATE="ESTIMATE" _LABEL="LABEL" _UPPER="UPPER" _LOWER="LOWER";
run;
ods pdf close;
