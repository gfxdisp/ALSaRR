# ALSaRR: Perceptual Model for Adaptive Local Shading and Refresh Rate

<img src="teaser.png"></img>

This repository contains the VRS motion quality dataset, source code for CaMoJAB metric and ALSaRR method described in our SIGGRAPH paper:

Akshay Jindal, Krzysztof Wolski, Karol Myszkowski and Rafał K. Mantiuk. 2021. Perceptual Model for Adaptive Local Shading and Refresh Rate. ACM Trans. Graph. (Proc. of SIGGRAPH Asia 2021), 40, 6, Article 280 (December 2021), 18 pages. https://doi.org/10.1145/3478513.3480514

The paper and videos can be found on the project web page: https://www.cl.cam.ac.uk/research/rainbow/projects/alsarr/

Repository structure (see individual file for usage and dependencies):

CaMoJAB:
    - camojab.m             -> content-aware motion quality metric
    - camojab_fits.mat      -> prefitted metric parameters
    - model_Dst.m           -> core logic of camojab metric
    - sample_texture.png
    - example.m             -> example usage of camojab metric
    - fits/train_camojab.m
    - fits/test_camojab.m
    - fits/ablation_camojab.m
    - paper_plots/plot_luminance_vs_persistence.m
    - paper_plots/plot_fps_vs_vel_vs_ppd.m
    - paper_plots/plot_fps_vs_ppd_vs_vel.m

ALSaRR:
    - adaptive_fps.m      -> (example) calculate optimal fps vs vel LUT
    - adaptive_vrs.m      -> (example) calculate optimal VRS map
    - vrs_dp.m            -> optimal solution to VRS knapsack problem
    - vrs_greedy.m        -> near optimal solution to VRS knapsack problem
    - polyfit_quality.m   -> (example) how to fit polynomial to camojab
    - sample_data/        -> sample data for above scripts

Experiment1_data:
    - vrs_exp_data.csv          -> Pairwise comparision data from Experiment 1
    - scale_pwc2jnd.m           -> Script to scale data to JND units
    - vrs_exp_data_jnd.mat      -> PWC data scaled to JND
    - vrs_exp_stimulus.mat      -> Textures used in the experiment
    - external_data/marrr_exp1_scaled.mat -> Data from [Denes et al. 2020] Experiment 1
    - external_data/marrr_exp3_scaled.mat -> Data from [Denes et al. 2020] Experiment 3
    - external_data/ITU.mat     -> Data from [BT Series. 2020] Fig. 7 and Fig. 20
    - external_data/mackins.mat -> Data from [Mackin et al. 2016] Fig. 2

Utils: collection of  helper functions 

