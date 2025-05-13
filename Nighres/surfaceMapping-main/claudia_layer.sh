#!/bin/sh

#for subject in *;

#do echo $subject;pwd
#cd $subject

 	rm layer_analysis/layer_mean.dat
	3dROIstats -mask path/to/layer.nii.gz -1DRformat -quiet -nzmean bold.nii.gz >>layer_analysis/layer_mean.dat
	rm layer_analysis/layer_mean_FINAL2.dat

	WRD=$(head -n 1 layer_analysis/layer_mean.dat|wc -w); for((i=2;i<=$WRD;i=i+2)); do awk '{print $'$i'}' layer_analysis/layer_mean2.dat| tr '\n' ' ';echo; done > layer_analysis/layer_mean_FINAL2.dat
 

#cd ../
#done