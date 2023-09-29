%macro mediansby;
	**OVERALL;
	proc summary noprint
		data=&dataset;
		var &var;
		where &var is not missing AND clinstatus ne 999;
		output out=outmedian n=n_all median=median_all q1=q1_all q3=q3_all;
	run;
	**BY CLINICAL STATUS;
	%do clinstatus=1 %to &&maxlevel;
	proc summary noprint
		data=&dataset;
		where clinstatus=&clinstatus AND &var is not missing;
		var &var;
		output out=outmedian_&clinstatus n=n_&clinstatus median=median_&clinstatus q1=q1_&clinstatus q3=q3_&clinstatus;
	run;
	%end;
	**P-value overall comparison;
	**P-values for post hoc comparisons;
/*	ods select dscf;*/
/*	ods output dscf=pairwise;*/
	proc npar1way kw dscf noprint
		data=&dataset;
		where &var is not missing AND clinstatus ne 999;
		class clinstatus;
		var &var;
		output out=mpval;
	run;
/*	data pairwise;*/
/*		set pairwise;*/
/*		if comparison in('1 vs. 2', '2 vs. 1') then pvalue_12=PValue;*/
/*			else if comparison in('1 vs. 3','3 vs. 1') then pvalue_13=Pvalue;*/
/*			else if comparison in('1 vs. 4','4 vs. 1') then pvalue_14=Pvalue;*/
/*			else if comparison in('2 vs. 3','3 vs. 2') then pvalue_23=Pvalue;*/
/*			else if comparison in('2 vs. 4','4 vs. 2') then pvalue_24=Pvalue; */
/*			else if comparison in('3 vs. 4','4 vs. 3') then pvalue_34=Pvalue;*/
/*			else if comparison in('1 vs. 5','5 vs. 1') then pvalue_15=Pvalue;*/
/*			else if comparison in('2 vs. 5','5 vs. 2') then pvalue_25=Pvalue;*/
/*			else if comparison in('3 vs. 5','5 vs. 3') then pvalue_35=Pvalue;*/
/*			else if comparison in('4 vs. 5','5 vs. 4') then pvalue_45=Pvalue;*/
/*	run;*/
	**Need to do separate wilcoxon noprint tests;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (1,2);
		class clinstatus;
		var &var;
		output out=mpval_12;
	run;
		proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (1,3);
		class clinstatus;
		var &var;
		output out=mpval_13;
	run;
		proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (1,4);
		class clinstatus;
		var &var;
		output out=mpval_14;
	run;
		proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (1,5);
		class clinstatus;
		var &var;
		output out=mpval_15;
	run;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (2,3);
		class clinstatus;
		var &var;
		output out=mpval_23;
	run;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (2,4);
		class clinstatus;
		var &var;
		output out=mpval_24;
	run;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (2,5);
		class clinstatus;
		var &var;
		output out=mpval_25;
	run;	
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (3,4);
		class clinstatus;
		var &var;
		output out=mpval_34;
	run;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (3,5);
		class clinstatus;
		var &var;
		output out=mpval_35;
	run;
	proc npar1way wilcoxon noprint
		data=&dataset;
		where clinstatus in (4,5);
		class clinstatus;
		var &var;
		output out=mpval_45;
	run;
	data mpval_12; set mpval_12; keep PT2_WIL; rename PT2_WIL=pvalue_12; run;
	data mpval_13; set mpval_13; keep PT2_WIL; rename PT2_WIL=pvalue_13; run;
	data mpval_14; set mpval_14; keep PT2_WIL; rename PT2_WIL=pvalue_14; run;
	data mpval_15; set mpval_15; keep PT2_WIL; rename PT2_WIL=pvalue_15; run;
	data mpval_23; set mpval_23; keep PT2_WIL; rename PT2_WIL=pvalue_23; run;
	data mpval_24; set mpval_24; keep PT2_WIL; rename PT2_WIL=pvalue_24; run;
	data mpval_25; set mpval_25; keep PT2_WIL; rename PT2_WIL=pvalue_25; run;
	data mpval_34; set mpval_34; keep PT2_WIL; rename PT2_WIL=pvalue_34; run;
	data mpval_35; set mpval_35; keep PT2_WIL; rename PT2_WIL=pvalue_35; run;
	data mpval_45; set mpval_45; keep PT2_WIL; rename PT2_WIL=pvalue_45; run;
	data pairwise; set mpval_:; run;
	proc summary
		data=pairwise;
		var pvalue_:;
		output out=mpval_post max=;
	run;
	data outmedians; 
		merge 	
			outmedian
			outmedian_1-outmedian_&&maxlevel
			mpval (keep=P_KW rename=(P_KW=pvalue))
			mpval_post
		; 
		length vartext $35. var $200. varlevel $8.;
		drop _type_ _freq_; 
		vartext="&var";
		var="(n="||strip(n_all)||")";
		varlevel="&varlevel";
		**CALCULATE IQR;
		array q1 (*) q1_all q1_1-q1_&&maxlevel;
		array q3 (*) q3_all q3_1-q3_&&maxlevel;
		array iqr (*) iqr_all iqr_1-iqr_&&maxlevel;
			do i=1 to dim(iqr);
				iqr(i)=q3(i)-q1(i);
			end;
		drop i;
	run;
	data outdatacontm; set outdatacontm outmedians; run;
%mend;
%macro LNmediansby;
	**OVERALL;
	proc summary noprint
		data=&dataset;
		var &var;
		where &var is not missing AND clinstatus ne 999;
		output out=outmedian n=n_all median=median_all q1=q1_all q3=q3_all;
	run;
	**BY CLINICAL STATUS;
	%do clinstatus=1 %to &&maxlevel;
	proc summary noprint
		data=&dataset;
		where clinstatus=&clinstatus AND &var is not missing;
		var &var;
		output out=outmedian_&clinstatus n=n_&clinstatus median=median_&clinstatus q1=q1_&clinstatus q3=q3_&clinstatus;
	run;
	%end;
	**LOG TRANSFORM FOR ALL COMPARISONS AND DO PARAMETRIC TEST;
	data ln;
		set &dataset;
		lnvar=log(&var);
	run;
	ods select Diff OverallANOVA;
	ods output OverallANOVA=LNmpval Diff=pairwise;
	proc glm
		data=ln;
		class clinstatus;
		where &var is not missing AND clinstatus ne 999;
		model lnvar=clinstatus;
		lsmeans clinstatus/pdiff;
	run;
	quit;
	data pairwise;
		set pairwise;
		if rowname=1 then do;
			pvalue_12=_2;
			pvalue_13=_3;
			pvalue_14=_4;
			pvalue_15=_5;
		end;
		if rowname=2 then do;
			pvalue_23=_3;
			pvalue_24=_4;
			pvalue_25=_5;
		end;
		if rowname=3 then do;
			pvalue_34=_4;
			pvalue_35=_5;
		end;
		if rowname=4 then do;
			pvalue_45=_5;
		end;
	run;
	proc summary
		data=pairwise;
		var pvalue_12 pvalue_13 pvalue_23 pvalue_14 pvalue_24 pvalue_34 pvalue_15 pvalue_25 pvalue_35 pvalue_45;
		output out=LNmpval_post max=;
	run;
	data LNmpval;
		set LNmpval;
		if source='Model';
		keep probF;
		rename probf=pvalue;
	run;
	data outmedians; 
		merge 	
			outmedian
			outmedian_1-outmedian_&&maxlevel
			LNmpval
			LNmpval_post
		; 
		length vartext $35. var $200. varlevel $8.;
		drop _type_ _freq_; 
		vartext="&var";
		var="(n="||strip(n_all)||")";
		varlevel="&varlevel";
		**CALCULATE IQR;
		array q1 (*) q1_all q1_1-q1_&&maxlevel;
		array q3 (*) q3_all q3_1-q3_&&maxlevel;
		array iqr (*) iqr_all iqr_1-iqr_&&maxlevel;
			do i=1 to dim(iqr);
				iqr(i)=q3(i)-q1(i);
			end;
		drop i;
	run;
	data outdatacontLNm; set outdatacontLNm outmedians; run;
%mend;
%macro freqsby;
	**OVERALL;
	proc freq noprint
		data=&dataset;
		where clinstatus ne 999 AND &var is not missing;
		table &var /out=outfreq;
	run;
	%do clinstatus=1 %to &&maxlevel;
	**BY CLINICAL STATUS;
	proc freq noprint
		data=&dataset;
		where clinstatus=&clinstatus AND &var is not missing;
		table &var/out=outfreq_&clinstatus;
	run;
	data outfreq_&clinstatus;
		set outfreq_&clinstatus;
		rename count=count_&clinstatus percent=percent_&clinstatus;
	run;
	%end;
	**GENERATE PVALUE;
	proc freq data=&dataset noprint; 
		table &var*clinstatus /chisq;
		output out=fpval chisq; 
	run;
	**GENERATE PVALUES FOR BETWEEN GROUP COMPARISONS;
	proc freq data=&dataset noprint; 
		where clinstatus in(1,2);
		table &var*clinstatus /chisq;
		output out=fpval_12 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(1,3);
		table &var*clinstatus /chisq;
		output out=fpval_13 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(2,3);
		table &var*clinstatus /chisq;
		output out=fpval_23 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(1,4);
		table &var*clinstatus /chisq;
		output out=fpval_14 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(2,4);
		table &var*clinstatus /chisq;
		output out=fpval_24 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(3,4);
		table &var*clinstatus /chisq;
		output out=fpval_34 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(1,5);
		table &var*clinstatus /chisq;
		output out=fpval_15 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(2,5);
		table &var*clinstatus /chisq;
		output out=fpval_25 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(3,5);
		table &var*clinstatus /chisq;
		output out=fpval_35 chisq; 
	run;
	proc freq data=&dataset noprint; 
		where clinstatus in(4,5);
		table &var*clinstatus /chisq;
		output out=fpval_45 chisq; 
	run;
	**Merge All;
	data outfreqs;
		length vartext $35. varlevel $8.;
		merge 	outfreq	(rename=(count=count_all percent=percent_all))
				outfreq_1-outfreq_&&maxlevel	
		;
		by &var;
		vartext="&var";
		if &var=99 then tagunknown=1;		**unknown;
		if &var=999 then tagunknown=2;		**missing;
		if &var=998 then tagunknown=3;		**unknown/missing;
		if &o=1 then tagunknown=_n_;
		var=vvalue(&var);
		varlevel="&varlevel";
		if &dropbin ne . then do;
			if &var=&dropbin then delete;
		end;
		drop &var;
	run;
	proc sort data=outfreqs; by tagunknown descending count_all; run;
	data fpval_12; set fpval_12; keep p_pchi; rename p_pchi=pvalue_12; run;
	data fpval_13; set fpval_13; keep p_pchi; rename p_pchi=pvalue_13; run;
	data fpval_14; set fpval_14; keep p_pchi; rename p_pchi=pvalue_14; run;
	data fpval_23; set fpval_23; keep p_pchi; rename p_pchi=pvalue_23; run;
	data fpval_24; set fpval_24; keep p_pchi; rename p_pchi=pvalue_24; run;
	data fpval_34; set fpval_34; keep p_pchi; rename p_pchi=pvalue_34; run;
	data fpval_15; set fpval_15; keep p_pchi; rename p_pchi=pvalue_15; run;
	data fpval_25; set fpval_25; keep p_pchi; rename p_pchi=pvalue_25; run;
	data fpval_35; set fpval_35; keep p_pchi; rename p_pchi=pvalue_35; run;
	data fpval_45; set fpval_45; keep p_pchi; rename p_pchi=pvalue_45; run;
	data outfreqs;
		merge outfreqs 	fpval (keep=P_PCHI rename=(P_PCHI=pvalue)) 
						fpval_:	;
	run;
	data outdatacat; set outdatacat outfreqs; run;
%mend;
%macro post;
	data outdata; 
		set outdatacat (in=a) /*outdatacont (in=b)*/ outdatacontm (in=c) outdatacontLNm (in=d); 
		if vartext='' & var='' then delete; 
		varlevel2=input(varlevel,best12.);
		drop varlevel; rename varlevel2=varlevel;
		if a then type=1;			/*freqs*/
			*else if b then type=2;	/*means*/
			else if c then type=3;	/*medians*/
			else if d then type=3;	/*medians + ln-transformed pvals*/
	run;
	proc sort data=outdata; by varlevel; run;
%mend;
