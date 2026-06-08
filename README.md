# llm_empathy_CiHB

This repository is for review process of the manuscript titled "Cross-Lingual Emotion Translation Using Multimodal Large Language Models: A Proof-of-Concept Study with Relevance to Psychopathological Profiling".

#1. Conda environment for analysis

```
conda env create --file llm_empathy_review_env.yaml
```

#2. Files that contain user responses : All items are de-identified.
- raw_data_emotions.csv : 547 English-speaking users' (Source cohort from manuscript) ratings of 18 emotions along core affect dimensions of valence, arousal, focus, dominance
- text_coreaffects.csv : 91 study participants' (Evaluation cohort) affective responses to 10 emotional texts per 18 emotions along core affect dimensions of valence, arousal, focus, dominance
- image_coreaffects.csv : 91 study participants' (Evaluation cohort) affective responses to 18 emotional images along core affect dimensions of valence, arousal, focus, dominance
- demographic_all2.csv : Demographic variables of Evaluation cohort participants
- factor_scores.csv : Pathology ratings of Evaluation cohort participants (PHQ-9, GAD-7, STAI-X1, ERQ, and PFQ-2)
- mug_korean_samples_affect.xslx : 76 Korean participants' (Validation cohort) ratings of 18 emotions along core affect dimensions of valence, arousal, focus, dominance
- df_pathology_magnitude.xlsx : Evaluation cohort's psychopathology and core affect magnitude values for Linear mixed effects model (R)
- df_pathology_redundancy.xlsx : Evaluation cohort's psychopathology and core affect redundancy values for Linear mixed effects model (R)

#3. Code script for replication of analyses
- llm_analysis_draft_chb.ipynb -> Affective transfer fidelity analysis in python
- lmm.R -> Psychopathology and affective magnitude/redundacny analysis in R
