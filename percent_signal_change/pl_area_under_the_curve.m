%% ........................................................................Tidy up
clear all
close all
clc
%% ......................................................................
load('D2_percent_signal_change_young.mat')

curve = [percent_signal_change_young(9,:);percent_signal_change_young(10,:);percent_signal_change_young(11,:);percent_signal_change_young(12,:);percent_signal_change_young(13,:);percent_signal_change_young(14,:);percent_signal_change_young(15,:)];

mean = nanmean(curve,2);

x = 1:7;
x = x';

D2_young = [];

for i_y = 1:size(curve,2)

    y = curve(:,i_y);

    area = trapz(x,y);

    D2_young = [D2_young; area];

end

D2_young_mean = trapz(x,mean);

%% ......................................................................
load('D2_percent_signal_change_old.mat')

curve = [percent_signal_change_old(9,:);percent_signal_change_old(10,:);percent_signal_change_old(11,:);percent_signal_change_old(12,:);percent_signal_change_old(13,:);percent_signal_change_old(14,:);percent_signal_change_old(15,:)];

mean = nanmean(curve,2);

x = 1:7;
x = x';

D2_old = [];

for i_y = 1:size(curve,2)

    y = curve(:,i_y);

    area = trapz(x,y);

    D2_old = [D2_old; area];

end

D2_old_mean = trapz(x,mean);

%% ......................................................................
load('D2_percent_signal_change_young.mat')

curve_1 = [percent_signal_change_young(7,:);percent_signal_change_young(8,:);percent_signal_change_young(9,:);percent_signal_change_young(10,:);percent_signal_change_young(11,:)];

mean_1 = nanmean(curve_1,2);

x = 1:5;
x = x';

area_under_the_curve_1 = [];

for i_y = 1:size(curve_1,2)

    y = curve_1(:,i_y);

    area = trapz(x,y);

    area_under_the_curve_1 = [area_under_the_curve_1; area];

end

mean_area_under_the_curve_1 = trapz(x,mean_1);

curve_2 = [percent_signal_change_young(14,:);percent_signal_change_young(15,:);percent_signal_change_young(16,:);percent_signal_change_young(17,:);percent_signal_change_young(18,:)];

mean_2 = nanmean(curve_2,2);

x = 1:5;
x = x';

area_under_the_curve_2 = [];

for i_y = 1:size(curve_2,2)

    y = curve_2(:,i_y);

    area = trapz(x,y);

    area_under_the_curve_2 = [area_under_the_curve_2; area];

end

mean_area_under_the_curve = trapz(x,mean_2);

D2_young = area_under_the_curve_1 + area_under_the_curve_2;
%% ......................................................................
load('D2+D3_percent_signal_change_old.mat')

curve_1 = [percent_signal_change_old(7,:);percent_signal_change_old(8,:);percent_signal_change_old(9,:);percent_signal_change_old(10,:);percent_signal_change_old(11,:)];

mean_1 = nanmean(curve_1,2);

x = 1:5;
x = x';

area_under_the_curve_1 = [];

for i_y = 1:size(curve_1,2)

    y = curve_1(:,i_y);

    area = trapz(x,y);

    area_under_the_curve_1 = [area_under_the_curve_1; area];

end

mean_area_under_the_curve_1 = trapz(x,mean_1);

curve_2 = [percent_signal_change_old(14,:);percent_signal_change_old(15,:);percent_signal_change_old(16,:);percent_signal_change_old(17,:);percent_signal_change_old(18,:)];

mean_2 = nanmean(curve_2,2);

x = 1:5;
x = x';

area_under_the_curve_2 = [];

for i_y = 1:size(curve_2,2)

    y = curve_2(:,i_y);

    area = trapz(x,y);

    area_under_the_curve_2 = [area_under_the_curve_2; area];

end

mean_area_under_the_curve = trapz(x,mean_2);

D2D3_old = area_under_the_curve_1 + area_under_the_curve_2;
