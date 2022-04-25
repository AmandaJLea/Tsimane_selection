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

We provide the code for the above analyses, but not the underlying genomic data. Genomic data collected from indigenous groups is subject to their sovereignty. This includes their determination of uses that benefit the communities and/or align with their priorities, that protect against the risk of identification of subjects or groups, and minimize uses with potentially stigmatizing interpretations. We are committed to open and reproducible science, but above all else we are committed to protecting the indigenous communities we work with and prioritizing their interests. The genomic data presented here were not collected under human subjects agreements or community agreements that explicitly discussed data sharing and secondary research; therefore, we cannot share genomic data at this time.

While we cannot make genomic data collected by the Tsimane Health and Life History Project (THLHP) available at this time, we are in the process of working with study communities to establish a formal process for data sharing and secondary research proposals. Information about this process will be posted to the THLHP website (https://tsimane.anth.ucsb.edu/index.html). 


Please contact me at amandalea7180@gmail.com with any questions.
