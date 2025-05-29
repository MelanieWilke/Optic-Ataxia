# Optic-Ataxia

## 1) code for Morel Atlas visualization
	
### Software used/required: 

```
[fList,pList] = matlab.codetools.requiredFilesAndProducts('ProcessMorelAtlas.m');
```
	
Name: 'MATLAB'                          Version: '9.13'
	
Name: 'Image Processing Toolbox'        Version: '11.6'
   	
### How to run:
- change to __\code_morel__ subfolder
- run __ProcessMorelAtlas.m__
	
to load and process the Morel Atlas data
	
which were prepared by Gunther Helms, PhD at the Department of Cognitive Neurology, University Medicine Goettingen, Germany
	
from the electronic version of the Morel Atlas (A. Krauth et al.)
	
https://www.research-collection.ethz.ch/handle/20.500.11850/669887
	
Persistent Link: https://doi.org/10.3929/ethz-b-000669887
	
made available under this licence: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International

- run __PlotMorelAtlas.m__ 
to generate a matlab figure of the Morel atlas on axial and coronal slices

[Settings in the beginning of the script control the details]


## 2) code for difference maps

under construction

### Software used/required
same as 1)

### How to run:
__currently disabled as code is without required behavioral data and lesion maps__
- change to __\code_morel__ subfolder
- run __plotVLSM_CombinedLesion.m__

general description:
 generate plots of voxelwise lesion overlaps and lesion subtraction maps
 overlayed on top of the MNI 0.5mm T1 template
 for all patients (general lesion overlap) and
 for patients with vs without a behavioral deficit.
 Images are decorated with overlayed outlines of selected thalamic ROIS of
 the Morel Atlas
 
 stores selected slices in axial and coronal orientation showing form left
 to right:
 1) lesion overlap for all patients with behavioral deficit
 2) lesion overlap for all patients without behavioral deficit
 3) lesion subtraction after transformation to % values (with - without
 deficits)
 4) lesion subtraction detail of thalamic region

 params.plotSlicesZ and params.plotSlicesY control which slices of the MNI template are shown
 
 iDeficit controls which behavioral variable is used and the cutoff for deficit
 
 The same code is used to generate a plot of lesion overlap for all patients (iDeficit = 5)
 
 Further settings in params

 iDrawingMode decides on whether all main lesions should be on one side (flipped), or whether thalamus masking is applied (only thalamic lesions are displayed)

 Also generates a summary image with multipe slices of the lesion subtraction details of the thalamic region variable "layout" controls the details of this plot 
 
 iVisual controls different visual styles for the plots

 All images are saved as matlab figures in \figures
