# README

This repository contains code for several analysis related to:

Natural selection of immune and metabolic genes associated with health in South American Amerindians

Amanda J. Lea, Angela Garcia, Jesusa Arevalo, Julien F. Ayroles, Kenneth Buetow, Steve W. Cole, Daniel Eid Rodriguez, Maguin Gutierrez, Heather M. Highland, Anne Justice, Thomas Kraft, Kari E. North, Jonathan Stieglitz, Hillard Kaplan, Benjamin C. Trumble, Michael D. Gurven

This manuscript uses an integrative dataset of genomic, transcriptomic, and immune and metabolic biomarker data for the Tsimane - an indigenous Amerindian population inhabiting the Bolivian lowlands - to identify loci putatively involved in local adaptation. 


The analyses that are described here are:

1) Processing Tsimane genotype data generated on the MEGA array for selection analyses (see genotype_processing.sh, PCRelate.R)
2) Calculating iHS, XP-EHH, and PBS from genotype data (see run_PBS_iHS_XP.sh) and combining these results to identify candidate regions for selection (see identify_candidate_regions.R)
4) Correlating SNPs within candidate regions with transcriptomic and phenotypic data (see RNAseq_processing.R, test_for_eQTL_LMM.R, test_for_associations_LMM.R)
5) Performing analyses of polygenic selection (see polygenic_selection_UKBB.R, polygenic_selection_VIPs.R)
6) Simulating genomic data generated from neutral demographic processes to compare to the observed data (see get_neutral_regions_for_demography.sh, demography.sh, compare_obs_sim_demography.R)


Data availability:

Individual-level data are stored in the Tsimane Health and Life History Project (THLHP) Data Repository, and are available through restricted access for ethical reasons. THLHP's highest priority is the safeguarding of human subjects and minimization of risk to study participants. The THLHP adheres to the ‘CARE Principles for Indigenous Data Governance’ (Collective Benefit, Authority to Control, Responsibility, and Ethics), which assure that the Tsimane and Moseten 1) have sovereignty over how data are shared, 2) are the primary gatekeepers determining ethical use, 3) are actively engaged in the data generation and 4) derive benefit from data generated and shared use whenever possible. The THLHP is also committed to the ‘FAIR Guiding Principles for scientific data management and stewardship’ (Findable, Accessible, Interoperable, Reusable). Requests for individual-level data should take the form of an application that minimally details the exact uses of the data and the research questions to be addressed, procedures that will be employed for data security and individual privacy, potential benefits to the study communities, and procedures for assessing and minimizing stigmatizing interpretations of the research results (see the following webpage for links to the data sharing policy and data request forms: https://tsimane.anth.ucsb.edu/data.html ). Requests for individual-level data will require institutional IRB approval (even if exempt) and will be reviewed by an Advisory Council composed of tribal leaders, tribal community members, Bolivian scientists, and the THLHP leadership. A similar structure exists for the Moseten data. The study authors and the Tsimane leadership are committed to open science and are available to assist interested investigators in preparing data access requests.

Please contact me at amandalea7180@gmail.com with any questions.
