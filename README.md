# EvaluatingEvasionStrategies
Supplementary Data and Code for PNAS manuscript of Jiao et al, 2023

This repository contains the data presented in the paper and the MATLAB code to process and analyze the data, as well as the MATLAB code to evaluate the probabilistic strategy models for evasion in zebrafish larvae.


Incomprehensive list of key files:


[dataAnalysis.m] data processing from pre-processed experiemntal measurements (not inucluded here) to the dataset we use for this study (evasionData.mat)

[NLL_\*.m] each computes the negative log likelihood (NLL) for a strategy model \n
[pX_\*.m] used by NLL*.m, each computes the PDF for a strategy model



\* includes DO, AP, CL. Note that OT and PL are incorporated into DO.

This repository is created and maintained by Yusheng Jiao, yusheng9559 at gmail.
