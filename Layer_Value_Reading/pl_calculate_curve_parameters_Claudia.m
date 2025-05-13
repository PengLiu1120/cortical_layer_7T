%% Derivative calculation
% This script reads in the h-file, binarizes the h-file, and multiplies the h-file with the T1-map, and reads out data

%% ........................................................................Tidy up
clear all
close all
clc

%% General Specifications
result_dir = '/Volumes/LayerPRF/LayerMapping/06_Derivative';

file_hand="allsubjects_T1_young_hand_short.mat";
% file_hand="allsubjects_T1_old_hand_short.mat";

n_layers = 21; % specifies number of layers used in final plots
start_layer = 1;
end_layer = 21;

n_fingers = 1;

%% load data table
load(file_hand);
table_allsubjects_hand = table_allsubjects_hand_short;
[num_rows, num_cols] = size(table_allsubjects_hand);

%% calculate mean layer profiles for all fingers and full finger map
mean_qT1_hand = nanmean(table_allsubjects_hand(:,1:n_layers),1);

%% get inflection points of qT1 1st derivative
% !!! only run with qT1 data, not with QSM data !!!

ind_inner= [];
ind_middle = [];
ind_outer = [];

maxima = zeros(n_fingers+1,n_layers);
minima = zeros(n_fingers+1,n_layers);
    
for j=1:n_fingers+2

    %calculate derivative vectors
    grad1=abs(gradient(mean_qT1_hand,1,2));
    grad2=gradient(grad1,1,2);
    % grad3=gradient(grad2,1,2);
    
    %find extrema of first derivative using matlab build-in functions
    local_minima = islocalmin(grad1,2);
    local_maxima = islocalmax(grad1,2);
    
    [r_max,c_max]=find(local_maxima);
    [r_min,c_min]=find(local_minima);
    
    maxima(j,c_max)=grad1(c_max);
    minima(j,c_min)=grad1(c_min);
    
    %find extrema using manual implementation
    %first look for exact zeros in 2nd derivative
    ind0_0 = find(grad2 == 0); 
    % then look for zero crossings between data points
    grad2_zerolocs = grad2(1:end-1) .* grad2(2:end);
    ind1_0 = find(grad2_zerolocs < 0);

    % zero locations of 2nd derivative
    ind = sort([ind0_0 ind1_0]);

    ind_inner = [ind_inner;3 ind(2)];
    ind_middle = [ind_middle;ind(2)+1 ind(3)];
    ind_outer = [ind_outer;ind(3)+1 21];
    
end

figure(1)
subplot(3,1,1)
plot(1:n_layers,mean_qT1_hand,'LineWidth',2,'Color','k')
ylim([1400 2400])
hold on;
subplot(3,1,2)
plot(1:n_layers,grad1,'LineWidth',2,'Color','k')
ylim([0 120])
hold on;
subplot(3,1,3)
plot(1:n_layers,grad2,'LineWidth',2,'Color','k')
ylim([-20 20])
hold on;