%let path=\\rc-fs.tch.harvard.edu\crc\Public\CRC General Information\Staff Folders\Milliren_Carly\Ozonoff - IMPACC\;
%include "&path\Reports\IMPACC Formats.sas";
%include "&path\Manuscript\IMPACC Clinical Manuscript DESC TABLE MACROS.sas";
options nofmterr;
data working;
	set traj_all;
	where visitday=0;
	**Core assay paper;
	if cohort_indicator=1;
run;

**************************;
**		TABLE 1			**;
**************************;
**COLUMNS - allow to change depending on definition of clinical status;
%let clinstatus=cluster5;
%let clinstatusname=Trajectory Groups;
%let maxlevel=5;					**levels in clinical status 3-5?;
%let startdata=working;				**dataset to start with; 
data runtable1;
	set &startdata;
	clinstatus=&clinstatus;
	if clinstatus in(999,.) then delete;
run;
%let dataset=runtable1;
%let tableout=outtable1;
proc freq data=&dataset; table &clinstatus /out=totcounts outcum; run;
	**Save n's to populate table headers;
data _null_; 
	set totcounts; 
	by &clinstatus;
	call symput(cats('groupn_',&clinstatus),strip(count));
	if last.&clinstatus then do;
		call symput('totaln',strip(cum_freq));
	end;
run;
**ROWS; 
data outdata; run; data outdatacat; run; data outdatacontm; run; data outdatacontLNm; run; /*data outdatacont; run;*/ 
%let tableout=table1					;
**DEMOGRAPHICS;
%let var=enrollage						;	%let varlevel=1; 										%mediansby;
%let var=sex							;	%let varlevel=2;		%let o=.;	%let dropbin=.;		%freqsby;
/*%let var=raceeth2						;	%let varlevel=4;		%let o=.;	%let dropbin=.;		%freqsby;*/
%let var=race							;	%let varlevel=4;		%let o=.;	%let dropbin=.;		%freqsby;
%let var=ethnicity						;	%let varlevel=5;		%let o=.;	%let dropbin=.;		%freqsby;
**COMORBIDITIES;	
%let var=comorb_htn 					;	%let varlevel=6;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_dm				;	%let varlevel=7;		%let o=.;	%let dropbin=0;		%freqsby;
/*%let var=comorb_anyresp				;	%let varlevel=8;		%let o=.;	%let dropbin=0;		%freqsby;*/
%let var=comorb_anyresp_NOasthma		;	%let varlevel=8.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_asthma			;	%let varlevel=8.2;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_cardiac 			;	%let varlevel=9;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_ckd 				;	%let varlevel=10;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_neoplasm			;	%let varlevel=11;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_neuro					;	%let varlevel=12;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_heme 			;	%let varlevel=13;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_rheum 			;	%let varlevel=14;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_isaric_liver			;	%let varlevel=16;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_hxtrans					;	%let varlevel=17;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_eversmkvape				;	%let varlevel=18;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_anysubstance			;	%let varlevel=18.5;		%let o=.;	%let dropbin=0;		%freqsby;
/*%let var=comorb_etoh 					;	%let varlevel=19;		%let o=.;	%let dropbin=0;		%freqsby;*/
/*%let var=comorb_thc_cbd 				;	%let varlevel=20;		%let o=.;	%let dropbin=0;		%freqsby;*/
/*%let var=comorb_drugs 					;	%let varlevel=21;		%let o=.;	%let dropbin=0;		%freqsby;*/
%let var=comorb_tb 						;	%let varlevel=23;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_hiv2 					;	%let varlevel=22;		%let o=.;	%let dropbin=0;		%freqsby;
/*%let var=comorb_isaric_other			;	%let varlevel=24;		%let o=.;	%let dropbin=0;		%freqsby;*/
%let var=comorb_endocrine				;	%let varlevel=24.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_psych					;	%let varlevel=24.2;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_thrombo					;	%let varlevel=24.3;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_gi						;	%let varlevel=24.4;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_obgyn					;	%let varlevel=24.5;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comorb_count					;	%let varlevel=24.6;										%mediansby;
%let var=comorb_countcat				;	%let varlevel=24.7;		%let o=1;	%let dropbin=.;		%freqsby;
%let var=bmi							;	%let varlevel=25;										%mediansby;
%let var=bmicat							;	%let varlevel=26;		%let o=1;	%let dropbin=.;		%freqsby;
%let var=comorb_obese_metabolic			;	%let varlevel=27;		%let o=.;	%let dropbin=0;		%freqsby;
**EXPOSURE;
%let var=exposure_hcw					;	%let varlevel=28; 		%let o=.;	%let dropbin=.; 	%freqsby;
%let var=exposure_facility				;	%let varlevel=29; 		%let o=.;	%let dropbin=.; 	%freqsby;
/*%let var=exposure_lab					;	%let varlevel=30; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
%let var=exposure_contact				;	%let varlevel=31; 		%let o=.;	%let dropbin=.; 	%freqsby;
**SYMPTOMS;		
%let var=timesymptom					;	%let varlevel=32;										%mediansby;
%let var=timesymptomcat					;	%let varlevel=33;		%let o=1;	%let dropbin=.;		%freqsby;
%let var=symptom_allresp 				;	%let varlevel=33.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_cough					;	%let varlevel=33.11;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_cough_prod				;	%let varlevel=33.12;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_haemoptysis			;	%let varlevel=33.13;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_sorethroat				;	%let varlevel=33.14;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_rhinorrhea				;	%let varlevel=33.15;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_wheeze					;	%let varlevel=33.16;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_chestpain				;	%let varlevel=33.17;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_dyspnea				;	%let varlevel=33.18;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_chestwall				;	%let varlevel=33.19;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_earpain				;	%let varlevel=33.191;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_conjunctivitis			;	%let varlevel=33.192;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_allsyst 				;	%let varlevel=33.2;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_fever 					;	%let varlevel=33.21;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_chills 				;	%let varlevel=33.22;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_shaking 				;	%let varlevel=33.23;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_fatigue 				;	%let varlevel=33.24;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_myalgia 				;	%let varlevel=33.25;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_lymphadenopathy		;	%let varlevel=33.26;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_allneuro 				;	%let varlevel=33.3;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_headache 				;	%let varlevel=33.31;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_confusion 				;	%let varlevel=33.32;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_seizure 				;	%let varlevel=33.33;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_syncope 				;	%let varlevel=33.34;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_anosmia				;	%let varlevel=33.35;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_allGI 					;	%let varlevel=33.4;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_abdpain 				;	%let varlevel=33.41;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_nausea 				;	%let varlevel=33.42;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_diarrhea				;	%let varlevel=33.43;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_allderm 				;	%let varlevel=33.5;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_skinrash 				;	%let varlevel=33.51;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_skinulcer				;	%let varlevel=33.52;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_allother				;	%let varlevel=33.6;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_arthralgia 			;	%let varlevel=33.61;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_edema 					;	%let varlevel=33.62;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_nonamb 				;	%let varlevel=33.63;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=symptom_bleeding				;	%let varlevel=33.64;	%let o=.;	%let dropbin=0;		%freqsby;
**ADMISSION CHARACTERISTICS;
%let var=hosp_admit_source3				;	%let varlevel=34;		%let o=.;	%let dropbin=.;		%freqsby;
%let var=hosp_admitlevel				;	%let varlevel=35; 		%let o=.;	%let dropbin=.; 	%freqsby;
%let var=hosp_adm_febrile				;	%let varlevel=36;		%let o=.;	%let dropbin=.;		%freqsby;
**CLINICAL STATUS AT BASELINE;
%let var=baseline_tx_icu				;	%let varlevel=37; 		%let o=.;	%let dropbin=.; 	%freqsby;
%let var=baseline_anyimg				;	%let varlevel=37.1;		%let o=.;	%let dropbin=.;		%freqsby;
%let var=baseline_img_infil				;	%let varlevel=37.2;		%let o=1;	%let dropbin=.;		%freqsby;
/*%let var=baseline_vented 				;	%let varlevel=38; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
/*%let var=baseline_tx_ecmo				;	%let varlevel=39; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
/*%let var=baseline_tx_nippv			;	%let varlevel=40; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
/*%let var=baseline_tx_o2therapy		;	%let varlevel=41; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
/*%let var=baseline_tx_o2mode			;	%let varlevel=42; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
/*%let var=baseline_tx_o2flow			;	%let varlevel=43; 		%let o=.;	%let dropbin=.; 	%freqsby;*/
%let var=baseline_highestresp			;	%let varlevel=44;		%let o=1;	%let dropbin=.;		%freqsby;
%let var=baseline_spfio2ratio_low		;	%let varlevel=46;										%mediansby;
%let var=baseline_spfio2ratio_lowcat2	;	%let varlevel=47;		%let o=1;	%let dropbin=.;		%freqsby;
%let var=baseline_sofa_score_new		;	%let varlevel=48;										%mediansby;
**BASELINE LABS;
%let var=baseline_lab_wbc				;	%let varlevel=49; 										%LNmediansby;
%let var=baseline_lab_lymph				;	%let varlevel=50; 										%LNmediansby;	
%let var=baseline_lymph_abn 			;	%let varlevel=50.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lymph_abn2 			;	%let varlevel=50.2;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_anc				;	%let varlevel=51; 										%LNmediansby;
%let var=baseline_anc_abn 				;	%let varlevel=51.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_anc_lymph_ratio		;	%let varlevel=51.5;										%LNmediansby;
%let var=baseline_lab_platelets			;	%let varlevel=52; 										%LNmediansby;
%let var=baseline_platelets_abn 		;	%let varlevel=52.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_platelets_abn2		;	%let varlevel=52.2;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_hgb				;	%let varlevel=53; 										%LNmediansby;
%let var=baseline_hgb_abn 				;	%let varlevel=53.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_pt				;	%let varlevel=54; 										%LNmediansby;
%let var=baseline_lab_inr				;	%let varlevel=55;										%LNmediansby;
%let var=baseline_lab_ast				;	%let varlevel=56; 										%LNmediansby;
%let var=baseline_lab_alt				;	%let varlevel=57; 										%LNmediansby;												
%let var=baseline_alt_abn				;	%let varlevel=57.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_bili				;	%let varlevel=58; 										%LNmediansby;
%let var=baseline_lab_cr				;	%let varlevel=59; 										%LNmediansby;
%let var=baseline_cr_abn 				;	%let varlevel=59.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_pct				;	%let varlevel=60; 										%LNmediansby;
%let var=baseline_lab_crp				;	%let varlevel=61; 										%LNmediansby;
%let var=baseline_crp_abn		 		;	%let varlevel=61.1;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=baseline_lab_ldh				;	%let varlevel=62; 										%LNmediansby;
%let var=baseline_lab_ddimer			;	%let varlevel=63; 										%LNmediansby;
%let var=baseline_ddimer_abn			;	%let varlevel=63.1;		%let o=.;	%let dropbin=0;		%freqsby; 
%let var=baseline_lab_ferritin			;	%let varlevel=64; 										%LNmediansby;
%let var=baseline_lab_trop				;	%let varlevel=65;										%LNmediansby;	
%let var=baseline_trop_abn				;	%let varlevel=65.1;		%let o=.;	%let dropbin=0;		%freqsby;
**LOS/Escalation of care;
%let var=los							;	%let varlevel=66;										%mediansby;
%let var=ever_ever_esc					;	%let varlevel=67;										%freqsby;
**MEDS;	
%let var=ever_conplas 					;	%let varnum=71;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_anyparasitic				;	%let varnum=72;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_remdes 					;	%let varnum=73;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antiIL6 					;	%let varnum=74;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antiIL1 					;	%let varnum=75;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_ccr5 						;	%let varnum=76;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_monoclonal 				;	%let varnum=77;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_jakinh 					;	%let varnum=78;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antiviral 				;	%let varnum=79;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_steroids 					;	%let varnum=80;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_vitamin 					;	%let varnum=81;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_azithro 					;	%let varnum=82;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_otherabx 					;	%let varnum=83;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_neuropsych 				;	%let varnum=84;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_nsaids 					;	%let varnum=85;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_aceinh 					;	%let varnum=86;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_arbs 						;	%let varnum=87;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_anticoag 					;	%let varnum=88;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_inhsteroid 				;	%let varnum=89;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antihist 					;	%let varnum=90;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antihyper 				;	%let varnum=91;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_hyperlipid 				;	%let varnum=92;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_immsupp 					;	%let varnum=93;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_pressors 					;	%let varnum=94;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_vaccine 					;	%let varnum=95;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antiarrhyth 				;	%let varnum=96;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_antifungal 				;	%let varnum=97;	%let o=.;	%let dropbin=0;		%freqsby;
/*%let var=ever_antiparasite 			;	%let varnum=98;	%let o=.;	%let dropbin=0;		%freqsby;*/
%let var=ever_pulmhyp 					;	%let varnum=99;	%let o=.;	%let dropbin=0;		%freqsby;
%let var=ever_othersupport				;	%let varnum=100;%let o=.;	%let dropbin=0;		%freqsby;
**COMPLICATIONS;	
%let var=anycomp						;	%let varlevel=101;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=ncomp							;	%let varlevel=102;										%mediansby;
%let var=comp_seizure 					;	%let varlevel=103;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_meningitis				;	%let varlevel=104;		%let o=.;	%let dropbin=0;		%freqsby; 
%let var=comp_stroke 					;	%let varlevel=105;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_afib					 	;	%let varlevel=106;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_ventricular_arr			;	%let varlevel=107;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_endo						;	%let varlevel=108;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_myocard 					;	%let varlevel=109;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_cardarrest 				;	%let varlevel=110;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_stemi 					;	%let varlevel=111;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_nstemi 					;	%let varlevel=112;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_chf 						;	%let varlevel=113;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_pneumonia 				;	%let varlevel=114;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_bronch 					;	%let varlevel=115;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_pneumo 					;	%let varlevel=116;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_pneumomed					;	%let varlevel=117;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_aat 						;	%let varlevel=118;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_avthromb					;	%let varlevel=119;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_shock 					;	%let varlevel=120;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_bac 						;	%let varlevel=121;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_coag 						;	%let varlevel=122;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_hypogly				 	;	%let varlevel=123;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_bleed 					;	%let varlevel=124;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_anemia 					;	%let varlevel=125;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_pancreas					;	%let varlevel=126;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_renal 					;	%let varlevel=127;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_liver 					;	%let varlevel=128;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_othlab 					;	%let varlevel=129;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_electrolyte 				;	%let varlevel=130;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_othinf 					;	%let varlevel=131;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_othpulm 					;	%let varlevel=132;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_unddisease 				;	%let varlevel=133;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_enceph 					;	%let varlevel=134;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_othneuro					;	%let varlevel=135;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_hypergly 					;	%let varlevel=136;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_psych 					;	%let varlevel=137;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_mis 						;	%let varlevel=138;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_gastroent 				;	%let varlevel=139;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_pulmvas					;	%let varlevel=140;		%let o=.;	%let dropbin=0;		%freqsby;
%let var=comp_other2					;	%let varlevel=141;		%let o=.;	%let dropbin=0;		%freqsby;
																										%post;	

options papersize=LETTER orientation=landscape topmargin="0.5in" bottommargin="0.5in" leftmargin="1in" rightmargin="1in" nodate center nonumber;
ods rtf file="&path\Manuscript\Output\IMPACC Manuscript Table CORE ASSAY COHORT by &clinstatusname &sysdate9..rtf";
ods escapechar="~";
title "Baseline Characteristics by &clinstatusname";

proc report
	data=outdata nowd split="#" SPANROWS
	style(report) = {font = ("calibri",10pt,normal)frame=hsides rules=all bordercolor=black}
	style(header) = {font = ("calibri",12pt,bold) foreground=black background = white}
	style(column) = {font = ("calibri",10pt,normal) foreground=black background = white};

	columns type ("	" (vartext var)) 
			count_all percent_all median_all iqr_all
			count_1 percent_1 median_1 iqr_1
			count_2 percent_2 median_2 iqr_2
			count_3 percent_3 median_3 iqr_3
			/**THREE GROUPS**/
/*			("n (%)" (all one two three))*/
/*			pvalue pvalue_12 pvalue_13 pvalue_23*/
/*			("Pairwise#group#comparisons" (comp_12 comp_13 comp_23))*/
			/**FOUR GROUPS**/
/*			count_4 percent_4 median_4 iqr_4*/
/*			("n (%)" (all one two three four))*/
/*			pvalue pvalue_12 pvalue_13 pvalue_23 pvalue_14 pvalue_24 pvalue_34*/
/*			("Pairwise#group#comparisons" (comp_12 comp_13 comp_14 comp_23 comp_24 comp_34))*/
			/**FIVE GROUPS**/
			count_4 percent_4 median_4 iqr_4 
			median_5 iqr_5 count_5 percent_5
			("n (%)" (all one two three four five))
			pvalue pvalue_12 pvalue_13 pvalue_23 pvalue_14 pvalue_24 pvalue_34 pvalue_15 pvalue_25 pvalue_35 pvalue_45
			("Pairwise#group#comparisons" (comp_12 comp_13 comp_14 comp_15 comp_23 comp_24 comp_25 comp_34 comp_35 comp_45))
	;
	define vartext 			/order order=data F=$vartext2f. " " style={font=("arial",10pt,bold)};
	define var 				/display " ";
	define type 			/analysis noprint;
	define count_all 		/analysis noprint;
	define percent_all		/analysis noprint;
	define median_all		/analysis noprint;
	define iqr_all			/analysis noprint;
	define all	 			/computed right "Overall#(N=&totaln)";
	define count_1 			/analysis noprint;
	define percent_1		/analysis noprint;
	define median_1			/analysis noprint;
	define iqr_1			/analysis noprint;
	define one	 			/computed right "GROUP 1#(n=&groupn_1)";
	define count_2 			/analysis noprint;
	define percent_2		/analysis noprint;
	define median_2			/analysis noprint;
	define iqr_2			/analysis noprint;
	define two	 			/computed right "GROUP 2#(n=&groupn_2)";
	define count_3 			/analysis noprint;
	define percent_3		/analysis noprint;
	define median_3			/analysis noprint;
	define iqr_3			/analysis noprint;
	define three	 		/computed right "GROUP 3#(n=&groupn_3)";
	define pvalue			/display analysis "Overall#p-value" F=pval.;
	define pvalue_12		/analysis noprint;
	define pvalue_13		/analysis noprint;
	define pvalue_23		/analysis noprint;
	define comp_12			/computed center "1 vs. 2";
	define comp_13			/computed center "1 vs. 3";
	define comp_23			/computed center "2 vs. 3";
	**ONLY NEEDED IF 4 GROUPS;
	define count_4 			/analysis noprint;
	define percent_4		/analysis noprint;
	define median_4			/analysis noprint;
	define iqr_4			/analysis noprint;
	define four	 			/computed right "GROUP 4#(n=&groupn_4)";
	define pvalue_14		/analysis noprint;
	define pvalue_24		/analysis noprint;
	define pvalue_34		/analysis noprint;
	define comp_14			/computed center "1 vs. 4";
	define comp_24			/computed center "2 vs. 4";
	define comp_34			/computed center "3 vs. 4";
/*	**ONLY NEEDED IF 5 GROUPS;*/
	define count_5 			/analysis noprint;
	define percent_5		/analysis noprint;
	define median_5			/analysis noprint;
	define iqr_5			/analysis noprint;
	define five	 			/computed right "GROUP 5#(n=&groupn_5)";
	define pvalue_15		/analysis noprint;
	define pvalue_25		/analysis noprint;
	define pvalue_35		/analysis noprint;
	define pvalue_45		/analysis noprint;
	define comp_15			/computed center "1 vs. 5";
	define comp_25			/computed center "2 vs. 5";
	define comp_35			/computed center "3 vs. 5";
	define comp_45			/computed center "4 vs. 5";

	**OVERALL;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute all /character length=30;
		**FREQS;
	if type.sum=1 & count_all.sum ne . then all=catt(put(count_all.sum,5.0),' (',compress(put(percent_all.sum,5.0)),'%)');
		else if type.sum=1 & count_all.sum=. then all="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_all.sum ne . then all=catt(put(median_all.sum,5.1),' (',compress(put(iqr_all.sum,5.1)),')');
		else if type.sum=3 & median_all.sum=. then all=" - ";
	endcomp;

	**GROUP 1;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute one /character length=30;
		**FREQS;
	if type.sum=1 & count_1.sum ne . then one=catt(put(count_1.sum,5.0),' (',compress(put(percent_1.sum,5.0)),'%)');
		else if type.sum=1 & count_1.sum=. then one="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_1.sum ne . then one=catt(put(median_1.sum,5.1),' (',compress(put(iqr_1.sum,5.1)),')');
		else if type.sum=3 & median_1.sum=. then one=" - ";
	endcomp;

	**GROUP 2;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute two /character length=30;
		**FREQS;
	if type.sum=1 & count_2.sum ne . then two=catt(put(count_2.sum,5.0),' (',compress(put(percent_2.sum,5.0)),'%)');
		else if type.sum=1 & count_2.sum=. then two="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_2.sum ne . then two=catt(put(median_2.sum,5.1),' (',compress(put(iqr_2.sum,5.1)),')');
		else if type.sum=3 & median_2.sum=. then two=" - ";
	endcomp;

	**GROUP 3;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute three /character length=30;
		**FREQS;
	if type.sum=1 & count_3.sum ne . then three=catt(put(count_3.sum,5.0),' (',compress(put(percent_3.sum,5.0)),'%)');
		else if type.sum=1 & count_3.sum=. then three="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_3.sum ne . then three=catt(put(median_3.sum,5.1),' (',compress(put(iqr_3.sum,5.1)),')');
		else if type.sum=3 & median_3.sum=. then three=" - ";
	endcomp;

	**PAIRWISE COMPARISONS;
		**THREE GROUPS;
	compute comp_12 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_12.sum<0.005 then comp_12="*";
		if .<pvalue_12.sum<0.001 then comp_12="**";
		if .<pvalue_12.sum<0.0001 then comp_12="***";
	end;
	endcomp;
	compute comp_13 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_13.sum<0.005 then comp_13="*";
		if .<pvalue_13.sum<0.001 then comp_13="**";
		if .<pvalue_13.sum<0.0001 then comp_13="***";
	end;
	endcomp;
	compute comp_23 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_23.sum<0.005 then comp_23="*";
		if .<pvalue_23.sum<0.001 then comp_23="**";
		if .<pvalue_23.sum<0.0001 then comp_23="***";
	end;
	endcomp;

	**ONLY NEEDED IF 4 GROUPS;
	**GROUP 4;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute four /character length=30;
		**FREQS;
	if type.sum=1 & count_4.sum ne . then four=catt(put(count_4.sum,5.0),' (',compress(put(percent_4.sum,5.0)),'%)');
		else if type.sum=1 & count_4.sum=. then four="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_4.sum ne . then four=catt(put(median_4.sum,5.1),' (',compress(put(iqr_4.sum,5.1)),')');
		else if type.sum=3 & median_4.sum=. then four=" - ";
	endcomp;
		**FOUR GROUPS;
	compute comp_14 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_14.sum<0.005 then comp_14="*";
		if .<pvalue_14.sum<0.001 then comp_14="**";
		if .<pvalue_14.sum<0.0001 then comp_14="***";
	end;
	endcomp;	
	compute comp_24 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_24.sum<0.005 then comp_24="*";
		if .<pvalue_24.sum<0.001 then comp_24="**";
		if .<pvalue_24.sum<0.0001 then comp_24="***";
	end;
	endcomp;
	compute comp_34 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_34.sum<0.005 then comp_34="*";
		if .<pvalue_34.sum<0.001 then comp_34="**";
		if .<pvalue_34.sum<0.0001 then comp_34="***";
	end;
	endcomp;

	*************************;
	**ONLY NEEDED IF 5 GROUPS;
	**GROUP 5;
	**N (%) or MEAN (SD) or MEDIAN (IQR);
	compute five /character length=30;
		**FREQS;
	if type.sum=1 & count_5.sum ne . then five=catt(put(count_5.sum,5.0),' (',compress(put(percent_5.sum,5.0)),'%)');
		else if type.sum=1 & count_5.sum=. then five="0 (0%)";
		**MEDIANS;
		else if type.sum=3 & median_5.sum ne . then five=catt(put(median_5.sum,5.1),' (',compress(put(iqr_5.sum,5.1)),')');
		else if type.sum=3 & median_5.sum=. then five=" - ";
	endcomp;
		**FIVE+ GROUPS;
	compute comp_15 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_15.sum<0.005 then comp_15="*";
		if .<pvalue_15.sum<0.001 then comp_15="**";
		if .<pvalue_15.sum<0.0001 then comp_15="***";
	end;
	endcomp;	
	compute comp_25 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_25.sum<0.005 then comp_25="*";
		if .<pvalue_25.sum<0.001 then comp_25="**";
		if .<pvalue_25.sum<0.0001 then comp_25="***";
	end;
	endcomp;
	compute comp_35 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_35.sum<0.005 then comp_35="*";
		if .<pvalue_35.sum<0.001 then comp_35="**";
		if .<pvalue_35.sum<0.0001 then comp_35="***";
	end;
	endcomp;
	compute comp_45 /character length=5;
	if .<pvalue.sum<0.05 then do;
		if .<pvalue_45.sum<0.005 then comp_45="*";
		if .<pvalue_45.sum<0.001 then comp_45="**";
		if .<pvalue_45.sum<0.0001 then comp_45="***";
	end;
	endcomp;

run;
footnote "Adjusted for multiple comparisons *p<0.05	**p<0.01	***p<0.001";
ods rtf close;

**adjust for multiple comparisons - Bonferroni;
/*	p<0.05=p<0.005*/
/*	p<0.01=p<0.001*/
/*	p<0.001=p<0.0001*/
