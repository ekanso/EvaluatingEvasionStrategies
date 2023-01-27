# EvaluatingEvasionStrategies
Supplementary Data and Code for PNAS manuscript of Jiao et al, 2023

This repository contains the data presented in the paper and the MATLAB code to process and analyze the data, as well as the MATLAB code to evaluate the probabilistic strategy models for evasion in zebrafish larvae.
## Data file
[evasionData.mat] provides all data used in this research project.  
other .mat files store the saved results from important/time-consuming evaluations.

## Incomprehensive list of key MATLAB scripts:
Sorry for the messy naming patterns...

[dataAnalysis.m] data processing from pre-processed experiemntal measurements (not inucluded here) to the dataset we use for this study (evasionData.mat)
[expDataHist.m] Histograms of the experimental data we used in the evaluations.  
[Fig_CurlToTurn.m] creates the figure that shows the mapping from curling angle to turning angle in zebrafish larvae.  
[Fig_heatmap.m] 2D histograms showing the joint distributions of each pair of kinematic variables.  
[DeterministicStrategyEval.m] evaluates the strategies with no sensorimotor noise, computes the KL divergence estimates between experimental data and model predictions.  
[DeterministicStrategy.m] Other analysis and plots about the evasion models with no noise.  
[probModelBootstrapping.m] bootstraps datasets for optimization of noise parameters in the probabilistic evasion models.  
[probModelEval.m] Visualize the evaluation results for the models with sensorimotor noise.  
[optimization_map/openloop/param_dist.m] files to illustrate the process of parameter optimization and boostrapping.  

[NLL_\*.m] each computes the negative log likelihood (NLL) for a strategy model.  
[pX_\*.m] used by NLL*.m, each computes the PDF for a strategy model.  
\* includes DO, AP, CL. Note that OT and PL are incorporated into DO.

[circular_linearFitting.m] customized linear regression method between two circular variables.  
[connectMatrix.m, getAddedMass.m, initialConfig.m, prescribedAngle.m, threelink_dynamics.m, simplifiedFishRun.m] scripts to run the three-link fish dynamics.  
[turningCompute.m] computes the turning angle by going around an ellipse trajectory on the shape space
[parameterSpace.m] visualizes the geometric colormap for analyzing the three-link fish model.



This repository is created and maintained by Yusheng Jiao, yusheng9559 at gmail.
