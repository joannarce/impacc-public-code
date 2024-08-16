%let filepath=\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\Manuscript;
%include "\\rc-fs.tch.harvard.edu\crc\public\crc general information\staff folders\milliren_carly\Ozonoff - IMPACC\Reports\IMPACC Formats.sas";
proc import
	datafile="&filepath\Convalescent Paper\Final Conv Clustering\IMPACC_ward3_6_10_data_06172022.xlsx"
	out=conv_cluster
	dbms=xlsx replace;
run;
proc import
	datafile="&filepath\Convalescent Paper\Final Conv Clustering\PROMIS_membership_06092022.xlsx"
	out=conv_traj
	dbms=xlsx replace;
	sheet="trajectory groups";
run;
data conv_traj2;
	set conv_traj;
	**Reverse code trajectories so higher is better;
	array traj (*) class_eq5d5l class_health class_cognitive class_dyspnea;
	array traj2 (*) class2_eq5d5l class2_health class2_cognitive class2_dyspnea;
		do i=1 to dim(traj);
			traj2(i)=4-traj(i);
		end;

	drop class_physical;
run;
proc freq data=conv_traj2; table class_eq5d5l class_health class_cognitive class_dyspnea; run;
%macro import;
proc import
	datafile="&filepath\Convalescent Paper\Final Conv Clustering\PROMIS_scores_05272022.xlsx"
	out=&out
	dbms=xlsx replace;
	sheet="&sheet";
run;
proc sort data=&out; by pid month; run;
%mend;
%let out=eq5d5l;		%let sheet=EQ-5D-5L; 		%import;
%let out=health;		%let sheet=Health;			%import;
%let out=psychosocial;	%let sheet=Impact;			%import;
%let out=mental;		%let sheet=Mental;			%import;
%let out=physical;		%let sheet=Physical;		%import;
%let out=cognitive;		%let sheet=Psychological;	%import;
%let out=dyspnea;		%let sheet=Dyspnea_new;		%import;
proc format;
	value clusterf 1='MIN' 2='PHY' 5='COG' 6='MLT';
run;
data all;
	merge eq5d5l health psychosocial mental physical cognitive dyspnea;
	by pid month;
	rename 
		PROMIS_Psycho=PROMIS_cognitive
		PROMIS_Impact=psychosocial_impact
	;
run;
proc sort data=conv_cluster; by pid; run;
data all_final; 
	merge conv_cluster all;
	by pid;
run;
proc freq data=conv_cluster; table cluster6; run;
proc summary
	data=all_final;
	class pid cluster6;
	var eq5d5l PROMIS_cognitive psychosocial_impact PROMIS_Dyspnea PROMIS_Mental health_score PROMIS_Physical;
	output out=summ mean=mean_eq5d5l mean_cognitive mean_psychosocial mean_dyspnea mean_mental mean_health mean_physical;
run;
proc format;
	value cluster6f
		1='MIN/cluster 1'
		2='MIN/cluster 3'
		3='MIN/cluster 4'
		4='PHY/cluster 2'
		5='COG/cluster 5'
		6='MLT/cluster 6'
	;
run;
data summ2;
	set summ;
	if _type_=3;
	drop _type_ _freq_;
	**Collapsed groupings;
	if cluster6 in(1,3,4) then cluster_final=1;
		else cluster_final=cluster6;
	format cluster_final clusterf.;

	**dummy vars for cluster;
	cluster_min=2; cluster_phy=2; cluster_cog=2; cluster_mlt=2;
	if cluster_final=1 then cluster_min=1;
	if cluster_final=2 then cluster_phy=1;
	if cluster_final=5 then cluster_cog=1;
	if cluster_final=6 then cluster_mlt=1;
	cluster2_1=cluster_min;
	cluster2_2=cluster_phy;
	cluster2_3=cluster_cog;
	cluster2_4=cluster_mlt;

		**dummy vars for original cluster;
	array clust (*) cluster_1-cluster_6;
		do i=1 to dim(clust);
			clust(i)=2;
			if cluster6=i then clust(i)=1;
		end;
	cluster6_new=cluster6;
	if cluster6=1 then cluster6_new=1;
		else if cluster6=3 then cluster6_new=2;
		else if cluster6=4 then cluster6_new=3;
		else if cluster6=2 then cluster6_new=4;
		else if cluster6=5 then cluster6_new=5;
		else if cluster6=6 then cluster6_new=6;

	format cluster6_new cluster6f.;
run;
proc sort data=conv_traj2; by pid; run;
data summ_final;
	merge summ2 conv_traj2; 
	by pid;
run;
proc means mean std median q1 q3 data=summ_final; var mean_physical mean_mental mean_psychosocial; run;
proc ttest
	data=summ_final;
	class cluster_mlt;
	var class2_health class2_eq5d5l mean_physical class2_dyspnea mean_mental mean_psychosocial class2_cognitive;
run;

%macro tests;
%do n=1 %to 6;
ods output TTests=tests_&n;
proc ttest
	data=summ_final;
	class cluster_&n;
	var class2_health class2_eq5d5l mean_physical class2_dyspnea mean_mental mean_psychosocial class2_cognitive;
run;
data tests_&n;
	set tests_&n;
	if method='Satterthwaite';
	keep variable tvalue probt;
	rename tvalue=tvalue_&n probt=probt_&n;
run;
%end;
%mend;
%tests;
data alltests;
	merge tests_1-tests_6;
run;
**TABLE 1AS - new table comparing 3 MIN clusters to justify collapsing them;
proc report
	data=alltests nowd split="#" SPANROWS;
		columns variable
				tvalue_1-tvalue_6 probt_1-probt_6
				min_1 min_3 min_4 
/*				phy_2 cog_5 mlt_6*/
		;
		define variable 	/order order=data;
		define tvalue_1 	/analysis noprint;
		define tvalue_2 	/analysis noprint;
		define tvalue_3 	/analysis noprint; 
		define tvalue_4 	/analysis noprint; 
		define tvalue_5 	/analysis noprint; 
		define tvalue_6 	/analysis noprint; 
		define probt_1 		/analysis noprint; 
		define probt_2 		/analysis noprint; 
		define probt_3 		/analysis noprint; 
		define probt_4 		/analysis noprint; 
		define probt_5 		/analysis noprint; 
		define probt_6 		/analysis noprint;
		define min_1 		/computed right "1/MIN";
		define min_3		/computed right "3/MIN";
		define min_4 		/computed right "4/MIN";
/*		define phy_2   		/computed right "2/PHY"; */
/*		define cog_5	 	/computed right "5/COG";*/
/*		define mlt_6 		/computed right "6/MLT";*/

		compute min_1 /character length=30;
			min_1=catt(put(tvalue_1.sum,8.2),' (',compress(put(probt_1.sum,pvalue.)),')');
		endcomp;
		compute min_3 /character length=30;
			min_3=catt(put(tvalue_3.sum,8.2),' (',compress(put(probt_3.sum,pvalue.)),')');
		endcomp;
		compute min_4 /character length=30;
			min_4=catt(put(tvalue_4.sum,8.2),' (',compress(put(probt_4.sum,pvalue.)),')');
		endcomp;
/*		compute phy_2 /character length=30;*/
/*			phy_2=catt(put(tvalue_2.sum,8.2),' (',compress(put(probt_2.sum,pvalue.)),')');*/
/*		endcomp;		*/
/*		compute cog_5 /character length=30;*/
/*			cog_5=catt(put(tvalue_5.sum,8.2),' (',compress(put(probt_5.sum,pvalue.)),')');*/
/*		endcomp;*/
/*		compute mlt_6 /character length=30;*/
/*			mlt_6=catt(put(tvalue_6.sum,8.2),' (',compress(put(probt_6.sum,pvalue.)),')');*/
/*		endcomp;*/
run;
**Pairwise comparisons between MIN for Table 1s;
%macro pair;
proc glm
	data=summ_final;
	class cluster6;
	where cluster6 in(1,3,4);
	model &var=cluster6;
	lsmeans cluster6 /pdiff /*adjust=tukey*/;
run;
quit;
%mend;
%let var=class2_health		;	%pair;
%let var=class2_eq5d5l		;	%pair;
%let var=mean_physical		;	%pair; 
%let var=class2_dyspnea 	;	%pair;
%let var=mean_mental 		;	%pair;
%let var=mean_psychosocial	;	%pair; 
%let var=class2_cognitive	;	%pair;





**TABLE 1BS (original Table 1S);
%macro tests;
%do n=1 %to 4;
ods output TTests=tests_&n;
proc ttest
	data=summ_final;
	class cluster2_&n;
	var class2_health class2_eq5d5l mean_physical class2_dyspnea mean_mental mean_psychosocial class2_cognitive;
run;
data tests_&n;
	set tests_&n;
	if method='Satterthwaite';
	keep variable tvalue probt;
	rename tvalue=tvalue_&n probt=probt_&n;
run;
%end;
%mend;
%tests;
data alltests;
	merge tests_1-tests_4;
run;

proc report
	data=alltests nowd split="#" SPANROWS;
		columns variable
				tvalue_1-tvalue_4 probt_1-probt_4
				min_1 phy_2 cog_5 mlt_6;
		;
		define variable 	/order order=data;
		define tvalue_1 	/analysis noprint;
		define tvalue_2 	/analysis noprint;
		define tvalue_3 	/analysis noprint; 
		define tvalue_4 	/analysis noprint; 
		define probt_1 		/analysis noprint; 
		define probt_2 		/analysis noprint; 
		define probt_3 		/analysis noprint; 
		define probt_4 		/analysis noprint; 
		define min_1 		/computed right "1/MIN";
		define phy_2   		/computed right "2/PHY"; 
		define cog_5	 	/computed right "5/COG";
		define mlt_6 		/computed right "6/MLT";

		compute min_1 /character length=30;
			min_1=catt(put(tvalue_1.sum,8.2),' (',compress(put(probt_1.sum,pvalue.)),')');
		endcomp;
		compute phy_2 /character length=30;
			phy_2=catt(put(tvalue_2.sum,8.2),' (',compress(put(probt_2.sum,pvalue.)),')');
		endcomp;		
		compute cog_5 /character length=30;
			cog_5=catt(put(tvalue_3.sum,8.2),' (',compress(put(probt_3.sum,pvalue.)),')');
		endcomp;
		compute mlt_6 /character length=30;
			mlt_6=catt(put(tvalue_4.sum,8.2),' (',compress(put(probt_4.sum,pvalue.)),')');
		endcomp;
run;
**Slight differences in the physical score but conclusion/result is the same;  
**Due to small difference in the data between mine and Shanshans versions;
**Didn't update the #s for physical score in the table for the revision to be consistent with the rest of the results - used
Shanshan's calculated t-statistics for Table 1BS;















**Plot individual outcomes;
%macro plot;
data plot;
	set summ2;
	var=&var;
run;
proc template;
define statgraph sgdesign;
dynamic _VAR _CLUSTER6_NEW;
begingraph / designwidth=706 designheight=545 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / 	xaxisopts=( type=linear display=(TICKS TICKVALUES LINE ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL )) 
						yaxisopts=( label=("&ylabel") labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL ));
         boxplot x=_CLUSTER6_NEW y=_VAR / name='box' display=(CAPS MEAN MEDIAN ) groupdisplay=Cluster outlineattrs=(color=CX000000 ) medianattrs=(color=CX000000 pattern=SOLID thickness=1 ) whiskerattrs=(color=CX000000 pattern=SOLID ) meanattrs=(color=CX000000 symbol=CIRCLEFILLED );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.plot template=sgdesign;
dynamic _VAR="'VAR'n" _CLUSTER6_NEW="'CLUSTER6_NEW'n";
run;
%mend;
%macro plot;
data plot;
	set summ2;
	var=&var;
run;
proc template;
define statgraph sgdesign;
dynamic _VAR _CLUSTER_FINAL;
begingraph / designwidth=706 designheight=545 border=false;
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay / 	xaxisopts=( type=discrete display=(TICKS TICKVALUES LINE ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL )) 
						yaxisopts=( label=("&ylabel") labelattrs=(color=CX000000 family='Arial' size=12 style=NORMAL weight=BOLD ) tickvalueattrs=(color=CX000000 family='Arial' size=10 style=NORMAL weight=NORMAL ));
         boxplot x=_CLUSTER_FINAL y=_VAR / name='box' display=(CAPS MEAN MEDIAN ) groupdisplay=Cluster outlineattrs=(color=CX000000 ) medianattrs=(color=CX000000 pattern=SOLID thickness=1 ) whiskerattrs=(color=CX000000 pattern=SOLID ) meanattrs=(color=CX000000 symbol=CIRCLEFILLED );
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.plot template=sgdesign;
dynamic _VAR="'VAR'n" _CLUSTER_FINAL="'CLUSTER_FINAL'n";
run;
%mend;
%let var=mean_health		;	%let ylabel=Health	; 				%plot;
%let var=mean_eq5d5l		;	%let ylabel=EQ-5D-5L;				%plot;	
%let var=mean_physical		;	%let ylabel=Physical;				%plot;
%let var=mean_dyspnea 		;	%let ylabel=Dyspnea;				%plot;
%let var=mean_mental 		;	%let ylabel=Mental;					%plot;
%let var=mean_psychosocial	;	%let ylabel=Psychosocial Impact;	%plot;
%let var=mean_cognitive		;	%let ylabel=Cognitive; 				%plot;
