# pupil_project
Code accompanying the paper: Chang, H.*, Tang, W.*, et al. (2024). Sleep micro-structure organizes memory replay. Nature.

+-------------------------------------------------------------------------+
|            Sleep micro-structure organizes memory replay                |
+-------------------------------------------------------------------------+
README.txt
Copyright (C) 2024, Wenbo Tang, version 1.0
All rights reserved.

DESCRIPTION OF THE CODE CONTAINED IN THE ARCHIVE: Tang_Nature_2024.tgz


BRIEF
=====

Code accompanying the paper: Chang, H.*, Tang, W.*, et al. (2024). Sleep micro-structure organizes memory replay. Nature.


GETTING STARTED
===============

Launch MATLAB and cd into the directory containing the code (e.g. '/pupil_project/Main_figure/').

Other files in the directory (with all sub-folders) needed in path:  
https://github.com/ayalab1/neurocode

Toolboxes required:
- Chronux (version 2.12; http://chronux.org/) 

- Uniform Manifold Approximation and Projection (UMAP) (version 4.1; https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

- structure_index (https://github.com/PridaLab/structure_index)

  For visualization:
- MatPlotLib (version 2.1.3; https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps)

- colormapformulae (version 1.0; https://www.mathworks.com/matlabcentral/fileexchange/66724-colormap-formulae) 

Data from animals included: HYC2, HYC3, PPP4, PPP7, PPP8, PPP10, PPP11, PPP12, PPP13, PPP14, PPP15, PPP16, PVR4, PPP24, PPP25

These codes were originally created in the MATLAB 2017a and 2023b. Source data produced from the scripts are provided along with the paper. 


FILES and FOLDERS
=================
  ./Fig1
  pupil_periodic_Fig1c.m : main script for calculating the periodicity of pupil-size dynamics in Fig. 1c
  sleepstate_umap_ALLstates_Fig1d.m : main script for generating example UMAP plots, demonstrating the sleep state-pupil relationship in Fig. 1d
  SleepState_LFPfeatures_ALLStates_preprocess_Fig1e.m: preprocessing script for gathering LFP features across all sleep states for calculating structure indices (SIs)
  SleepState_LFPfeatures_NREM_preprocess_Fig1e.m: preprocessing script for gathering LFP features within NREM for calculating structure indices
  SleepState_rawLFPfeatures_ALLStates_StructureIndex_Fig1e.m: main script for calculating SI using LFP features across all sleep states shown in Fig. 1e
  SleepState_rawLFPfeatures_NREM_StructureIndex_Fig1e.m: main script for calculating SI using LFP features within NREM shown in Fig. 1e
  RippleLFP_pupil_umap_example_Fig1g.m:  main script for generating example UMAP plots demonstrating the ripple-pupil relationship in Fig. 1g
  RippleLFP_pupil_prerprocessing_Fig1h.m: preprocessing script gathering ripple LFP features during NREM sleep for calculating structure indices (SIs)
  RippleLFP_pupil_rawStructureIndx_Fig1h.m: main script for calculating SI using ripple LFP features within NREM shown in Fig. 1h

  ./Fig2
  Replay_triggered_pupil_acrossanimal_matchMUA_Fig2b.m  : main script for calculating replay probabilities across sextiles of pupil size during NREM shown in Fig. 2b. Firing-rate matched distribution shown in EDFig. 5d is also calculated. See also .../src/ReplayDecoding/ for replay decoding demonstration
  AllSWR_properties_triggered_pupil_forGLM_Fig2ce.m: preprocessing script for gathering ripple and replay properties in relation to pupil size for GLMs and measuring replay percentage. 
  Pupil_replayprc_stats_Fig2c.m	: main script for calculating replay percentages across sextiles of pupil sizes during NREM shown in Fig. 2c
  glmgain_ripple_predict_pupil_Fig2e.m: main script for predicting pupil size using ripple and replay properties shown in Fig. 2e

  ./Fig3
  cal_Cheeseboard_pathLength_Fig3b.m : main script for calculating path length in the Cheeseboard task shown in Fig. 3. Other behavioral measures are provided in Source Data files.

  ./Fig4 
  Reactivation_FamiliarNovel_triggered_pupil_acrossanimal_Fig4e.m : main script for calculating reactivation strength in novel versus familiar environments across pupil sizes, shown in Fig. 4e. See also .../src/ReactivationStrength/ for demonstration of calculating RS. 
  Preexist_Reactivation_triggered_pupil_acrossanimal_Fig4g.m : main script for calculating reactivation strength for rigid (pre-existing) vs. plastic assemblies across pupil sizes, shown in Fig. 4g
  RankOrder_pupil_acrossanimal_Fig4h : main script for rank order correlation between PRE and POST sleep cross pupil sizes, shown in Fig. 4h

  ./Fig5
  MUAINT_triggered_pupil_sleep_acrossanimal_Fig5b.m : main script for measuring the relationship between pupil size and PYR and INT firing rate shown in Fig. 5b and EDFig. 10
  CCGtransmission_pupilRippletile_acrossanimals_Fig5d.m	: main script for gathering spike transmission probability of monosynaptic PYR-INT pairs across pupil sizes, shown in Fig. 5d
  CCGtransmission_pupilRippletile.m : function for calculating spike transmission probability of monosynaptic PYR-INT pairs across pupil sizes
  OptoRippleamp_SpontRipple_triggered_pupil_acrossanimal_Fig5f.m: main script for calculating amplitude of spontaneous ripples across pupil sizes shown in Fig. 5f. 
  OptoRippleamp_triggered_pupil_acrossanimal_Fig5f.m: main script for calculating amplitude of optogenetically induced ripples across pupil sizes shown in Fig. 5f. 

  ./src : 
  subfolder containing demonstrations of key analysis steps, including detecting SWR-associated high-synchrony events (HSEs), calculating reactivation strength, and replay decoding. 

  ./utilities : 
  subfolder containing all the helper functions.

  ./toolbox : 
  subfolder containing all the toolboxes used. Please refer to the original work from the developers.


CITING OUR WORK
===============

If you find the code useful, please cite the code source and the paper:
    Chang, H.*, Tang, W.*, Wulf, A. M., Nyasulu, T., Wolf, M. E., Fernandez-Ruiz, A., & Oliva, A.,(2024). Sleep micro-structure organizes memory replay. Nature.


CONTACT
=======
Bug reports, comments and questions are appreciated.
Please write to: 
	Wenbo Tang <wenbo.tang07@gmail.com>
