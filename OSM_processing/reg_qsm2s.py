import ants

subj = "ID"
subj1 = "ID"

out_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/qsm/"
anat_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/output/PrepareWBData/exp-0000/exp-0000-AAA/Intensity_Bounds/"
mov_dir = "/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/qsm/use_002_def_msdi2_230424_020924/" #"/media/juli/WIP-B2-SSD/onehander/data/" + subj + "/processed/output/reorient_function/exp-0000/exp-0000-A/reorient/"

fixed = subj1 + "_7-t1_mp2rage-foci_sag_p2_0_transform_calc_bound.nii.gz"
moving = subj1 + "_sensemap3_qsm_INTEGRAL_2_MSDI_multiplied_clone_transform.nii.gz"

reg = "reg_qsm2s.txt"

interpolation_method = "nearestNeighbor"

# Apply transformation
warpedImage = ants.apply_transforms(
	fixed = ants.image_read(anat_dir + fixed),
	moving = ants.image_read(mov_dir + moving),
	transformlist = (out_dir + reg),
	interpolator = interpolation_method,
	verbose = True,
	)

# Write file to disk
ants.image_write(warpedImage, out_dir + subj1 + "_sensemap3_qsm_INTEGRAL_2_MSDI_multiplied_clone_transform_reg.nii.gz") 
