import ants

subj = "ID" 
subj1 = "ID" 

out_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/" + subj1 + "_localisers/"
anat_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/output/PrepareWBData/exp-0000/exp-0000-AAA/Intensity_Bounds/"
#ROI_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/masks/"
loc_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/output/reorient_function/exp-0000/"
mov_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/" + subj1 + "_localisers/" #"/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/output/reorient_function/exp-0000/exp-0000-A/reorient/"

fixed = subj1 + "_7-t1_mp2rage-foci_sag_p2_0_transform_calc_bound.nii.gz"
moving = "data_Tmean_brain_clone_transform.nii.gz" #"data_Tmean_brain_resampled_clone_transform.nii.gz"

init_reg = 'f2s_init.txt'
#ROI = 'S1.left_right_final.nii.gz'

localizers_dirs = ["exp-0000-B/reorient/","exp-0000-C/reorient/","exp-0000-D/reorient/","exp-0000-E/reorient/","exp-0000-F/reorient/","exp-0000-G/reorient/","exp-0000-H/reorient/","exp-0000-I/reorient/","exp-0000-J/reorient/"]

#localizers = ["spmT_0002_left_hand_clone_transform.nii.gz","spmT_0002_left_hand_0_01_k3_peakc_bin_clone_transform.nii.gz","spmT_0003_right_hand_clone_transform.nii.gz","spmT_0004_left_foot_clone_transform.nii.gz","spmT_0004_left_foot_0_01_k3_peakc_bin_clone_transform.nii.gz","spmT_0005_right_foot_clone_transform.nii.gz","spmT_0005_right_foot_0_01_k0_peakc_bin_clone_transform.nii.gz","spmT_0006_tongue_clone_transform.nii.gz","spmT_0006_tongue_0_01_k3_peakc_bin_clone_transform.nii.gz"]
localizers = ["spmT_0002_left_hand_clone_transform.nii.gz","spmT_0003_right_hand_clone_transform.nii.gz","spmT_0003_right_hand_0_01_k3_peakc_bin_clone_transform.nii.gz","spmT_0004_left_foot_clone_transform.nii.gz","spmT_0004_left_foot_0_01_k3_peakc_bin_clone_transform.nii.gz","spmT_0005_right_foot_clone_transform.nii.gz","spmT_0005_right_foot_0_01_k3_peakc_bin_clone_transform.nii.gz","spmT_0006_tongue_clone_transform.nii.gz","spmT_0006_tongue_0_01_k3_peakc_bin_clone_transform.nii.gz"]

interpolation_method = "nearestNeighbor"

registeredImage = ants.registration(
	fixed = ants.image_read(anat_dir + fixed),
	moving = ants.image_read(mov_dir + moving),
	initial_transform = (out_dir + init_reg),
	#mask = ants.image_read(ROI_dir + ROI),
	#mask_all_stages = True,
	type_of_transform = 'SyNBoldAff',
	#aff_random_sampling_rate = 1,
	# aff_smoothing_sigmas = (3,2,1,0),
	#grad_step = 0.1,
	reg_iterations = (200, 100, 50 , 50),
	verbose = 1
	#outprefix = out_dir + "f2s_init" + '_final.txt',
	)

# Apply transformation
warpedImage = ants.apply_transforms(
	fixed = ants.image_read(anat_dir + fixed),
	moving = ants.image_read(mov_dir + moving),
	transformlist = registeredImage['fwdtransforms'],
	interpolator = interpolation_method,
	verbose = True,
	)

# Write file to disk
ants.image_write(warpedImage, out_dir + "data_Tmean_brain_clone_transform_reg" + '.nii.gz') 

# Apply transformation
#warpedImage = ants.apply_transforms(
#	fixed = ants.image_read(anat_dir + fixed),
#	moving = ants.image_read(mov_dir + moving),
#	transformlist = registeredImage['invtransforms'],
#	interpolator = interpolation_method,
#	verbose = True,
#	)

# Write file to disk
#ants.image_write(warpedImage, out_dir + "data_Tmean_brain_clone_transform_inv_reg" + '.nii.gz') 

i = 0

for localizer in localizers:

	# Apply transformation
	warpedImage = ants.apply_transforms(
		fixed = ants.image_read(anat_dir + fixed),
		moving = ants.image_read(loc_dir + localizers_dirs[i] + localizer),
		transformlist = registeredImage['fwdtransforms'],
		interpolator = interpolation_method,
		verbose = True,
		)

	# Write file to disk
	ants.image_write(warpedImage, out_dir + "analyses/" + localizer + '_reg.nii.gz') 

	i = i + 1

#ants.write_transform(registeredImage['fwdtransforms'][0], out_dir + "f2s_init_final.txt")

