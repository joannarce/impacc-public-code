This is a readme document for file ‘impacc_conv_symptoms.csv’.

Symptoms reported at baseline (on hospital admission) or during the convalescent period are categorized into three mutually exclusive groups:
- Cardiopulmonary: dyspnea, cough;
- Systemic: fatigue, myalgia, chills, fever;
- Other: nausea/vomiting (GI); anosmia, headache (neurologic), sore throat, red eyes/conjunctivitis (upper respiratory).

This file has 9 columns:1) participant_id.  Numeric identifier which is unique to each study participant.  Format is 0XX-YZZZ where X=site identifier (01-15), Y=site location (for sites which recruit participants across multiple locations), ZZZ=three digit sequential identifier (001-199).  Healthy controls use format CX where X takes values 6-10.2) convsx_cardio_1plus.  Indicates symptoms associated with cardiopulmonary system reported at 1+ timepoints during convalescence.  Values:
	0 = cardiopulmonary symptom not reported at 1+ convalescent timepoints;
	1 = cardiopulmonary symptom reported at 1+ convalescent timepoints;
	NA = data unavailable.
3) convsx_systemic_1plus.  Indicates symptoms associated with systemic illness reported at 1+ timepoints during convalescence.  Values:
	0 = systemic symptom not reported at 1+ convalescent timepoints;
	1 = systemic symptom reported at 1+ convalescent timepoints;
	NA = data unavailable.
4) convsx_other_1plus.  Indicates symptoms associated with other system (GI, neuro, or upper respiratory) reported at 1+ timepoints during convalescence.  Values:
	0 = other symptom not reported at 1+ convalescent timepoints;
	1 = other symptom reported at 1+ convalescent timepoints;
	NA = data unavailable.
5) convsx_1plus_count.  A count variable summing the variables for each of three system groups (cardiopulmonary, systemic, or other) indicating a symptom reported at 1+ convalescent timepoints.  This variable takes integer values from 0-3 or NA when data not available.
6) convsx_cardio_2plus.  Indicates symptoms associated with cardiopulmonary system reported at 2+ timepoints, including baseline.  Values:
	0 = cardiopulmonary symptom not reported at 2+ convalescent timepoints;
	1 = cardiopulmonary symptom reported at 2+ convalescent timepoints.
7) convsx_systemic_2plus.  Indicates symptoms associated with systemic illness reported at 2+ timepoints, including baseline.  Values:
	0 = systemic symptom not reported at 2+ convalescent timepoints;
	1 = systemic symptom reported at 2+ convalescent timepoints.
8) convsx_other_2plus.  Indicates symptoms associated with other system (GI, neuro, or upper respiratory) reported at 2+ timepoints, including baseline.  Values:
	0 = other symptom not reported at 2+ convalescent timepoints;
	1 = other symptom reported at 2+ convalescent timepoints.
9) convsx_2plus_count.  A count variable summing the variables for each of three system groups (cardiopulmonary, systemic, or other) indicating a symptom reported at 2+ timepoints, including baseline.  This variable takes integer values from 0-3.
