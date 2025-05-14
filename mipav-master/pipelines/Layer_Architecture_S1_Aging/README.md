This directory includes the structural MR image preprocessing and layering pipeline used in Liu, Doehler et al. 2023 (https://doi.org/10.1101/2023.12.01.567841) which is based on MIPAV (v7.3.0), JIST (v2.0-2013Nov06) and cbstools (3.0).

Installation instructions of these tools are provied in our Wiki.

This pipeline was particularly developed for processing of high resolution MP2RAGE data (0.5 mm). It consists of 5 major parts and 4 additional mapping pipelines:

01_Registration: co-registration of MP2RAGE whole brain and part brain images.\
02_MergingSkullStripping: merging of MP2RAGE whole brain and part brain images, skull and dura removal\
03_Segmentation: segmenting the brain into grey matter, white matter, cerebrospinal fluid\
04_PrepareMappingData: applying preprocessing steps to the original unchanged data (reduces interpolation steps to a minimum)\
05_Mapping_WholeBrain: cortical layering and cortical depth-dependent value sampling, mapping of sampled whole brain data onto the cortical surface (output in vtk format)\

06_MappingLeftS1: mapping of sampled data restricted to the left S1 region (can be applied to other ROIs)\
07_Mapping_QSM_WholeBrain: mapping of signed QSM whole brain data onto the cortical surface (can be applied to other layer-specific MRI contrasts)\
08_Mapping_function_without_layering: mapping of layer-unspecific functional (BOLD fMRI) data (1.0 mm resolution) onto the cortical surface (can be applied to other layer-unspecific MRI contrasts)\
09_Extract_Cortical_Thickness: layer-specific surface-based cortical thickness extraction using a quatitative T1 contrast

In addition further example pipelines to map (layer-unspecific) functional data onto cortical surfaces are provided:\
Mapping_function_ecm\
Mapping_function_hand_face_motor_maps\
Mapping_function_prf_without_layering

Example data to test the pipeline can be provided on request.