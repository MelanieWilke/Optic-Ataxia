# Optic-Ataxia
1) code for difference maps

under construction

2) code for Morel Atlas visualization
	
	Software used/required:
		[fList,pList] = matlab.codetools.requiredFilesAndProducts('VisualizeMorelAtlas.m');
	      Name: 'MATLAB'
          Version: '9.13'
          Name: 'Image Processing Toolbox'
          Version: '11.6'
    
	
	How to run:
	- change to \code subfolder
	- run ProcessMorelAtlas.m
		to load and process the Morel Atlas data
		
		which were prepared by Gunther Helms, PhD 
		at the Department of Cognitive Neurology, University Medicine Goettingen, Germany
		
		from the electronic version of the Morel Atlas (A. Krauth et al.)
		https://www.research-collection.ethz.ch/handle/20.500.11850/669887
		Persistent Link
		https://doi.org/10.3929/ethz-b-000669887
		made available under this licence:
		Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
		
	- run PlotMorelAtlas.m 
		to generate a matlab figure of the Morel atlas in axial and coronal slices
		[Settings in the beginning of the script control the details]
