%% Aperture vectorise
% This script takes the old aperture file and vectorise to fit in SamSrf ver 8.0 onwards
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 17 Oct 2022
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

%%
ApsFile = 'aps_verbars';

Scaling = [1 1];

VectoriseApertures(ApsFile, Scaling)

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));