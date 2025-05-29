% plotVLSM_CombinedLesion.m
% generate plots of voxelwise lesion overlaps and lesion subtraction maps
% overlayed on top of the MNI 0.5mm T1 template
% for all patients (general lesion overlap) and
% for patients with vs without a behavioral deficit.
% Images are decorated with overlayed outlines of selected thalamic ROIS of
% the Morel Atlas
% 
% stores selected slices in axial and coronal orientation showing form left
% to right:
% 1) lesion overlap for all patients with behavioral deficit
% 2) lesion overlap for all patients without behavioral deficit
% 3) lesion subtraction after transformation to % values (with - without
% deficits)
% 4) lesion subtraction detail of thalamic region
%
% params.plotSlicesZ and params.plotSlicesY control which slices  
% of the MNI template are shown
% iDeficit controls which behavioral variable is used and the cutoff for
% deficit
% The same code is used to generate a plot of lesion overlap for all
% patients (iDeficit = 5)
% further settings in params
%
% iDrawingMode decides on whether all main lesions should be on one side
% (flipped), or whether thalamus masking is applied (only thalamic lesions
% are displayed)
%
% Also generates a summary image with multipe slices of the lesion
% subtraction details of the thalamic region
% variable "layout" controls the details of this plot
% iVisual controls different visual styles for the plots
%
% all images are saved as matlab figures in \figures


% 20250512: change color bar for combined plots to positive only and
% include absolute numbers on color scales

clear;
close all;

LEFT = 1;
RIGHT = 2;

% check file was started correctly from \code_lesionmap
% generate paths to dataFolder and figFolder
thisFolder = pwd;
lastFilesepPosPlus1 = regexp(thisFolder,'code_lesionmap$');
if ~isempty(lastFilesepPosPlus1)
    codeFolder = thisFolder;
    lastFilesepPos = lastFilesepPosPlus1 - 1; % pointer was at 'c', now is at filesep
    % generate name of dataFolder
    dataFolder = [codeFolder(1:lastFilesepPos), 'data'];
    % generate name of figFolder
    figFolder = [codeFolder(1:lastFilesepPos), 'figures'];
else
    error('please start ProcessMorelAtlas from ...\Optic-Ataxia\code_lesionmap')
end



%% outer loop: choose iDeficit: possible values: 1:5
for iDeficit = 1:5
    switch iDeficit
        case 1 % sensory deficit, coded as 0 or -1 (deficit)
            params.criterionVariable ='CL_C_Sensory';
            params.criterionThreshold = -0.5;
        case 2 % motor deficit, coded as 0 or -1 (deficit)
            params.criterionVariable ='CL_C_Motor';
            params.criterionThreshold = -0.5;
        case 3 % contralesional hand and space optic ataxia
            params.criterionVariable ='AT_Error_CH_CS_Ataxia_HCRandom_z';
            params.criterionThreshold = -2;
        case 4 % smarties test deficit, grasping
            params.criterionVariable ='SM_Error_CH_C_Fov_z';
            params.criterionThreshold = -2;
        case 5 % all lesions
            params.criterionVariable ='CL_C_Sensory'; %with criterionTreshold == inf: this is all lesions
            params.criterionThreshold = inf;
    end

    % inner loop: choose iDrawingMode
    % paper uses modes 1 and 3: flipped or unflipped restricted to thalamus
    for iDrawingMode = 2:2:4
        switch iDrawingMode
            case 1
                params.doFlipping = true;
                params.doThalamusMasking = true;
            case 2
                params.doFlipping = true;
                params.doThalamusMasking = false;
            case 3
                params.doFlipping = false;
                params.doThalamusMasking = true;
            case 4
                params.doFlipping = false;
                params.doThalamusMasking = false;
        end


        % parameter for drawing
        % plotting of individual slices
        params.replaceColorMap = true;

        % for the single slice plots restricted to the slides shown in the
        % paper. For more complete slices use commented out part
        % slices are MNI 152 0.5mm slice numbers
        % which are later transformed to MNI coordinates
        params.plotSlicesZ = [154, 160]; %[10:20:130, 132:2:180, 190:20:350];
        params.plotSlicesY = [214]; %[10:20:170, 182:2:230, 240:20:430];
        
        % threshold for difference maps
        params.cutoffDifferenceImage = 0.21;
        params.colorSaturationAt = 1;
        if iDeficit == 5 % settings for all lesions
            params.cutoffDifferenceImage = 0.03; % show any lesion
            params.colorSaturationAt = 0.75;
        end
        
        % colorScheme was changed to make images more colorblind-friendly
        params.colorScheme = 'colorblind';%'jet'
        
        % Outlines of Morel ROIs
        params.ROIColors = 255; % color for general regions 255 (white), 0 (black)
        params.RegionDetailsMonochrome = true;
        params.showOnlyRegionsDetails = true;
        % ROIs to show on whole brain plots
        params.showRegions = {'MNI_L_Lateral_group'; 'MNI_R_Lateral_group'; 'MNI_L_Medial_group'; 'MNI_R_Medial_group';'MNI_L_Anterior_group'; 'MNI_R_Anterior_group'; 'MNI_L_Posterior_group'; 'MNI_R_Posterior_group'};
        % ROIs to show on detail plots
        params.showRegionsDetail = sort({'MNI_L_Pu.'; 'MNI_R_Pu.';'MNI_L_VA.'; 'MNI_R_VA.'; 'MNI_L_VL.'; 'MNI_R_VL.'; 'MNI_L_VM.'; 'MNI_R_VM.'; 'MNI_L_VPL.'; 'MNI_R_VPL.'; 'MNI_L_VPM.'; 'MNI_R_VPM.';'MNI_L_MD.'; 'MNI_R_MD.';'MNI_L_CM.'; 'MNI_R_CM.';'MNI_L_CL.'; 'MNI_R_CL.'});
        params.useNonOverlapForDrawing = true;
        nROIsDetail = numel(params.showRegionsDetail)/2;
        
        % Morel region gray-values for use in plots
        textcol = {};
        textstr = {};
        for iROI = 1:numel(params.showRegionsDetail)
            %hsvcol = [(rem(iROI-1,nROIsDetail)+1)/(nROIsDetail+1),0.75,1];
            %col = hsv2rgb(hsvcol); %[rem(rem(iROI,nROIs)/nROIs,1), rem(0.33 +rem(iROI,nROIs)/nROIs,1), rem(0.66 +rem(iROI,nROIs)/nROIs,1)] ;
            textcol{rem(iROI-1,nROIsDetail)+1} = 0.75*(rem(iROI-1,nROIsDetail)+1)/numel(params.showRegionsDetail)*[1, 1, 1];
            temp = regexp(params.showRegionsDetail {iROI}, '_','split');
            textstr{rem(iROI-1,nROIsDetail)+1} = temp{end}(1:end-1);
        end

        % study table with balanced HC for ataxia task and z transformation
        studyTableFileName = [dataFolder, filesep, 'ThalamusData_Ataxia_HCRandom_13_Feb_2023_15_41_01_CS_short.xlsx'];
        studyTable = readtable(studyTableFileName);

        % MNI template 0.5mm
        MNIImageFileName = [dataFolder, filesep, 'MNI152_T1_0.5mm.nii.gz'];
        mni = niftiread(MNIImageFileName);
        mniHeader = niftiinfo(MNIImageFileName);

        % Morel Atlas folder
        morelAtlasFolder = [dataFolder, filesep, 'morel-atlas_2024_09_16'];
        morelAtlasWeightsFolder = [dataFolder, filesep, 'morel-atlas_2024_09_16_weights'];
        
        % Lesions folder
        lesionFolder = [dataFolder, filesep, 'LesionImagesBrainRenamed'];
        lesionFolderFlipped = [dataFolder, filesep, 'LesionImagesBrainRenamedFlippedFlirt'];
        %lesionFolderForVLSM = 'D:\ThalamusGenerateDataTable\data\LesionImagesSameSideForVLSM';
        lesionFolderMasked = [dataFolder, filesep, 'LesionImagesBrainRenamedMasked'];
        lesionFolderFlippedMasked = [dataFolder, filesep, 'LesionImagesBrainRenamedFlippedFlirtMasked'];

        % reduce studyTable to patients ('Nxxx')
        subjects = studyTable.Subject;
        iPatients = cellfun(@(x) startsWith(x, 'N'), subjects);
        iHC = cellfun(@(x) startsWith(x, 'K'), subjects);
        HCTable = studyTable(iHC,:);
        studyTable = studyTable(iPatients,:);
        studyTable.Properties.RowNames = subjects(iPatients);
        subjects = studyTable.Subject;
        lesionSide = studyTable.LES_ThalamusSideFromLesion_HCRandom;

        % calculate new variables
        studyTable.CL_C_Sensory = zeros(size(studyTable.CL_C_Sens_Lower));
        studyTable.CL_C_Sensory(studyTable.CL_C_Sens_Lower==1 | studyTable.CL_C_Sens_Upper==1) = -1;
        studyTable.CL_C_Motor = zeros(size(studyTable.CL_C_Motor_Lower));
        studyTable.CL_C_Motor(studyTable.CL_C_Motor_Lower==1 | studyTable.CL_C_Motor_Upper==1) = -1;
        studyTable.AT_Error_Min_Congruent_Ataxia_HCRandom_z = min(studyTable.AT_Error_CH_CS_Ataxia_HCRandom_z, studyTable.AT_Error_IH_IS_Ataxia_HCRandom_z);


        % find lesion filenames corresponding to the patients in the study
        lesionFilenames = ExtractAndFilterFiles(lesionFolder, 'fnirt-05mm-wholebrain.nii.gz');
        lesionFilenamesFlipped = ExtractAndFilterFiles(lesionFolderFlipped, 'fnirt-05mm-wholebrain.nii.gz');

        iLesionFilenames = cellfun(@find, cellfun(@(x) cellfun(@(y) contains(y,x), lesionFilenames), subjects, UniformOutput=false));
        lesionFilenames = lesionFilenames(iLesionFilenames);
        iLesionFilenamesFlipped = cellfun(@find, cellfun(@(x) cellfun(@(y) contains(y,x), lesionFilenamesFlipped), subjects, UniformOutput=false));
        lesionFilenamesFlipped = lesionFilenamesFlipped(iLesionFilenamesFlipped);

        % morel regions
        morelRegionNames = ExtractAndFilterFiles(morelAtlasFolder,'.nii.gz');
        morelRegionWeightsNames = ExtractAndFilterFiles(morelAtlasWeightsFolder,'.nii.gz');
        morelHierarchyLevels = readtable([dataFolder, filesep, 'Morel_hierarchical_2025_05.xlsx']);

        ROINames = {};
        ROI = {};
        if params.useNonOverlapForDrawing
            % new 20240925: use weight == 1 to find exclusive ROIs
            for iROI = 1: numel(params.showRegions)
                ROINames{iROI,1} = morelRegionWeightsNames{contains(morelRegionWeightsNames, regexprep(params.showRegions{iROI},'\.','_'))};
                ROI{iROI,1} = niftiread(ROINames{iROI,1});
                ROI{iROI,1} = ROI{iROI,1}==1;
            end
        else
            for iROI = 1: numel(params.showRegions)
                ROINames{iROI,1} = morelRegionNames{contains(morelRegionNames, params.showRegions{iROI})};
                ROI{iROI,1} = niftiread(ROINames{iROI,1});
            end
        end

        ROINamesDetail = {};
        ROIDetail = {};
        if params.useNonOverlapForDrawing
            for iROI = 1: numel(params.showRegionsDetail)
                ROINamesDetail{iROI,1} = morelRegionWeightsNames{contains(morelRegionWeightsNames, regexprep(params.showRegionsDetail{iROI},'\.','_'))};
                ROIDetail{iROI,1} = niftiread(ROINamesDetail{iROI,1});
                ROIDetail{iROI,1} = ROIDetail{iROI,1}==1;
            end
        else
            for iROI = 1: numel(params.showRegionsDetail)
                ROINamesDetail{iROI,1} = morelRegionNames{contains(morelRegionNames, params.showRegionsDetail{iROI})};
                ROIDetail{iROI,1} = niftiread(ROINamesDetail{iROI,1});
            end
        end


        % morelMask
        morelMask = false(size(mni));
        for iROI = 1: numel(morelRegionNames)
            morelMask = morelMask | logical(niftiread(morelRegionNames{iROI}));
        end

        
        % thalamic (Morel) lesions
        nLesionsThal = zeros(size(mni), 'single');
        LesionsThal = zeros(size(mni),'uint32');
        for i = 1: numel(lesionFilenames)
            if lesionSide(i) == 1
                temp = niftiread(lesionFilenamesFlipped{i}).*morelMask;
            else
                temp = niftiread(lesionFilenames{i}).*morelMask;
            end
            nLesionsThal = nLesionsThal + temp;
            LesionsThal = LesionsThal + uint32(2^(i-1)) * uint32(temp);

        end
        nVoxelsAffectedInPatientsThal = arrayfun(@(x) sum(nLesionsThal(:)==x), 1:max(nLesionsThal(:)));
        uniqueLesionPatternsThal = unique(LesionsThal(:));
        nUniqueLesionPatternsThal = numel(uniqueLesionPatternsThal);
        nPatientsLesionedFromUniqueLesionPatternsThal = arrayfun(@(x) numel(strfind(dec2bin(x),'1')), uniqueLesionPatternsThal);
        nVoxelsLesionedFromUniqueLesionPatternsThal = arrayfun(@(x) sum(LesionsThal(:)==x), uniqueLesionPatternsThal);
        % unique logical vectors
        % convert integers to binary char array
        uniqueLesionPatternsBinThal = dec2bin(uniqueLesionPatternsThal);
        % convert to logical array
        % rows are subjects (after transpose)
        % (first subject is in row 1 after reversing uniqueLesionPatternsBinThal)
        % columns are patterns of lesions
        LesionVectorsThal = (uniqueLesionPatternsBinThal(:,end:-1:1) == '1')';

        nLesMin = 8;
        nLesMax = numel(lesionFilenames)-nLesMin + 1;
        iValidPatternsThal = find(nPatientsLesionedFromUniqueLesionPatternsThal >= nLesMin & nPatientsLesionedFromUniqueLesionPatternsThal <= nLesMax);
        nPerm = 10000;
        permutation = randperm(size(studyTable,1))
        behaviour = studyTable.(params.criterionVariable);
        behaviour = behaviour(permutation)

        BMLesionPatternThal = nan(size(uniqueLesionPatternsThal,1), nPerm + 1);
        for p = 1:nPerm + 1
            behaviour = behaviour(permutation);
            if p > 1
                permutation = randperm(size(studyTable,1));
                behaviour = studyTable.(params.criterionVariable);
            end
            for i = 1:numel(iValidPatternsThal)
                lesionStatus = LesionVectorsThal(:,iValidPatternsThal(i));
                lesionedGroupBehaviour = behaviour(lesionStatus);
                nonLesionedGroupBehaviour = behaviour(~lesionStatus);
                lesionedGroupBehaviourMatrix = repmat(lesionedGroupBehaviour,1,numel(nonLesionedGroupBehaviour));
                nonLesionedGroupBehaviourMatrix = repmat(nonLesionedGroupBehaviour', numel(lesionedGroupBehaviour),1);
                score = (0.5*(lesionedGroupBehaviourMatrix == nonLesionedGroupBehaviourMatrix) + (lesionedGroupBehaviourMatrix > nonLesionedGroupBehaviourMatrix))-0.5;
                meanScore = mean(score,'all');
                BMLesionPatternThal(iValidPatternsThal(i),p) = meanScore;
            end
        end





        % lesions
        iBelowCriterion = studyTable.(params.criterionVariable) <= params.criterionThreshold & cellfun(@(x) startsWith(x, 'N'),subjects);
        subjectsBelowCriterion = subjects(iBelowCriterion);
        iAboveCriterion = studyTable.(params.criterionVariable) > params.criterionThreshold & cellfun(@(x) startsWith(x, 'N'),subjects);
        subjectsAboveCriterion = subjects(iAboveCriterion);
        lesionSideBelowCriterion = lesionSide(iBelowCriterion);
        lesionSideAboveCriterion = lesionSide(iAboveCriterion);
        FileNamesBelowThreshold = cellfun(@(x) [lesionFolder, filesep,  x '_TwoThirdxff-to-MNI-fnirt-05mm-wholebrain.nii.gz'], subjectsBelowCriterion, 'UniformOutput', false);
        FileNamesAboveThreshold = cellfun(@(x) [lesionFolder, filesep,  x '_TwoThirdxff-to-MNI-fnirt-05mm-wholebrain.nii.gz'], subjectsAboveCriterion, 'UniformOutput', false);
        FileNamesBelowThresholdFlipped = cellfun(@(x) [lesionFolderFlipped, filesep,  x '_flipFlirtedBin_TwoThirdxff-to-MNI-fnirt-05mm-wholebrain.nii.gz'], subjectsBelowCriterion, 'UniformOutput', false);
        FileNamesAboveThresholdFlipped = cellfun(@(x) [lesionFolderFlipped, filesep,  x '_flipFlirtedBin_TwoThirdxff-to-MNI-fnirt-05mm-wholebrain.nii.gz'], subjectsAboveCriterion, 'UniformOutput', false);
        fprintf('%i subjects below criterion:\n', numel(FileNamesBelowThreshold))

        % load masks below criterion (deficit)
        gbc = single(zeros(size(mni)));
        for i = 1: numel(FileNamesBelowThreshold)
            FileNamesBelowThreshold{i}
            if params.doFlipping && lesionSideBelowCriterion(i) == 1
                % flip left lesions
                gbc = gbc + niftiread(FileNamesBelowThresholdFlipped{i});
            else
                gbc = gbc + niftiread(FileNamesBelowThreshold{i});
            end
        end
        if params.doThalamusMasking
            gbc = gbc.*morelMask;
        end

        % load masks above criterion (no deficit)
        fprintf('%i subjects above criterion:\n', numel(FileNamesAboveThreshold))
        gac = zeros(size(mni));
        for i = 1: numel(FileNamesAboveThreshold)
            FileNamesAboveThreshold{i}
            if params.doFlipping && lesionSideAboveCriterion(i) == 1
                % flip left lesions
                gac = gac + niftiread(FileNamesAboveThresholdFlipped{i});
            else
                gac = gac + niftiread(FileNamesAboveThreshold{i});
            end
        end
        if params.doThalamusMasking
            gac = gac.*morelMask;
        end

        % prepare colormap for lesion overlap and lesion difference
        colorMap = colormap(jet);
        % 202504: new colorMap to make it more red-green-blind-friendly
        if strcmp(params.colorScheme, 'colorblind')
            colorMap2 = zeros(size(colorMap));
            red = colorMap2(:,1);
            red(129:256)= linspace(0.5,1,128);
            red(128-7:128) = linspace(0,0.5,8);
            green = colorMap2(:,2);
            green(1:128) = linspace(1,0,128);
            green(129:256)= linspace(0,1,128);
            blue = colorMap2(:,3);
            blue(1:128) = linspace(0.5,0.8,128);
            blue(129:129+7) = linspace(0.8,0,8);
            colorMap = [red,green,blue];
        end
        
        colorMapScale = linspace(-1,1,256)';

        % prepare new colorMapNew with the extreme values
        colorMapReplace = [repmat(colorMap(1,:),128,1); repmat(colorMap(end,:),128,1)];
        % replace the inner part with a stretched version of the original map
        iColorMapReplace = colorMapScale >= -params.colorSaturationAt & colorMapScale <= params.colorSaturationAt;
        numColorMapReplace = sum(iColorMapReplace);
        colorMapScaleReplace = linspace(-1, 1, numColorMapReplace)';
        for iChannel = 1:3
            colorMapReplace(iColorMapReplace,iChannel) = interp1(colorMapScale, colorMap(:,iChannel), colorMapScaleReplace);
        end
        % replace colorMap With colorMapNew
        if params.replaceColorMap
            colorMap = colorMapReplace;
        end


        % relative numbers
        % gbc group below criterion
        % use monochrome MNI to initialize all three color channels
        redGBC = mni;
        greenGBC = mni;
        blueGBC =  mni;
        for i = 1:max(unique(gbc))
            round(128.5 + 127.5*i/numel(FileNamesBelowThreshold))
            redGBC(gbc==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesBelowThreshold)),1));
            greenGBC(gbc==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesBelowThreshold)),2));
            blueGBC(gbc==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesBelowThreshold)),3));
        end


        % gac group above criterion
        redGAC = mni;
        greenGAC = mni;
        blueGAC =  mni;
        for i = 1:max(unique(gac))
            round(128.5 + 127.5*i/numel(FileNamesAboveThreshold))
            redGAC(gac==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesAboveThreshold)),1));
            greenGAC(gac==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesAboveThreshold)),2));
            blueGAC(gac==i) = round(256 * colorMap(round(128.5 + 127.5*i/numel(FileNamesAboveThreshold)),3));
        end

        % gbcminusgac difference image
        if numel(FileNamesAboveThreshold)>0
            gbcminusgac = double(gbc)/numel(FileNamesBelowThreshold)-double(gac)/numel(FileNamesAboveThreshold);
        else
            gbcminusgac = double(gbc)/numel(FileNamesBelowThreshold)-double(gac);
        end
        redGBCminusGAC = mni;
        greenGBCminusGAC = mni;
        blueGBCminusGAC =  mni;
        uniquegbcminusgac = unique(gbcminusgac);
        for i = 1:numel(uniquegbcminusgac)
            % apply threshold
            if abs(uniquegbcminusgac(i)) >= params.cutoffDifferenceImage
                redGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),1));
                greenGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),2));
                blueGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),3));
            end
        end




        %% plotSlices Axial
        for useSlice =  params.plotSlicesZ
            MNI_Z = mniHeader.Transform.T(4,3) + (useSlice-0.5)* mniHeader.Transform.T(3,3);
            showImageGBC = uint8(zeros(size(mni,1), size(mni,2), 3));
            showImageGBC(:,:,1) = squeeze(redGBC(:,:,useSlice));
            showImageGBC(:,:,2) = squeeze(greenGBC(:,:,useSlice));
            showImageGBC(:,:,3) = squeeze(blueGBC(:,:,useSlice));
            showImageGBC(size(mni,1)+50, size(mni,2),3) = 0;
            % colorbar
            colorbarYPosGBC = size(mni,1)+10:size(mni,1)+20;
            colorbarXPosGBC = size(mni,2)/2:size(mni,2)-20;
            showImageGBC((colorbarYPosGBC(1)-1):(colorbarYPosGBC(end)+1),(colorbarXPosGBC(1)-1):(colorbarXPosGBC(end)+1),:) = 255;

            % map to [0:1]
            colorbarXPercentGBC = 100*(colorbarXPosGBC-min(colorbarXPosGBC))/(max(colorbarXPosGBC)-min(colorbarXPosGBC));
            colorbarXRangeGBC = 128.5 + 127.5 * colorbarXPercentGBC/100;
            colorbarXColorsGBC = round(256*(interp1(1:256,colorMap,colorbarXRangeGBC)));
            for iYPos = 1: numel(colorbarYPosGBC)
                showImageGBC(colorbarYPosGBC(iYPos),colorbarXPosGBC,:) = colorbarXColorsGBC;
            end
            % percent scale
            colorbarXPercentLabelsGBC = 0:10:100;
            colorbarXPercentLabelsPosGBC = colorbarXPosGBC(arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXPercentLabelsGBC));
            showImageGBC((colorbarYPosGBC(end)+1):(colorbarYPosGBC(end)+5),colorbarXPercentLabelsPosGBC,:) = 255;

            % subjects scale
            stepsizeGBC = ceil(numel(subjectsBelowCriterion)/10);
            colorbarXSubjectLabelsGBC = 1:stepsizeGBC:numel(subjectsBelowCriterion);
            colorbarXSubjectLabelsGBCfine = 1:numel(subjectsBelowCriterion)
            colorbarXSubjectPercentLabelsGBC = 100*colorbarXSubjectLabelsGBC/numel(subjectsBelowCriterion);
            colorbarXSubjectPercentLabelsGBCfine = 100*colorbarXSubjectLabelsGBCfine/numel(subjectsBelowCriterion);
            temp = arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXSubjectPercentLabelsGBC, 'UniformOutput', false);
            tempfine = arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXSubjectPercentLabelsGBCfine, 'UniformOutput', false);
            colorbarXSubjectPercentLabelsPosGBC = colorbarXPosGBC(cellfun(@(x) x(1), temp));
            colorbarXSubjectPercentLabelsPosGBCfine = colorbarXPosGBC(cellfun(@(x) x(1), tempfine));

            showImageGBC((colorbarYPosGBC(1)-5):(colorbarYPosGBC(1)-1),colorbarXSubjectPercentLabelsPosGBC,:) = 255;
            showImageGBC(colorbarYPosGBC,(colorbarXPosGBC(1)):(colorbarXSubjectPercentLabelsPosGBC(1)-1),:) = 0;
            showImageGBC(colorbarYPosGBC,(colorbarXSubjectPercentLabelsPosGBCfine(max(unique(gbc)))+1):(colorbarXPosGBC(end)),:) = 0;



            showImageGAC = uint8(zeros(size(mni,1), size(mni,2), 3));
            showImageGAC(:,:,1) = squeeze(redGAC(:,:,useSlice));
            showImageGAC(:,:,2) = squeeze(greenGAC(:,:,useSlice));
            showImageGAC(:,:,3) = squeeze(blueGAC(:,:,useSlice));
            showImageGAC(size(mni,1)+50, size(mni,2),3) = 0;
            if numel(subjectsAboveCriterion)>0
                % colorbar
                colorbarYPosGAC = size(mni,1)+10:size(mni,1)+20;
                colorbarXPosGAC = size(mni,2)/2:size(mni,2)-20;
                showImageGAC((colorbarYPosGAC(1)-1):(colorbarYPosGAC(end)+1),(colorbarXPosGAC(1)-1):(colorbarXPosGAC(end)+1),:) = 255;

                % map to [0:1]
                colorbarXPercentGAC = 100*(colorbarXPosGAC-min(colorbarXPosGAC))/(max(colorbarXPosGAC)-min(colorbarXPosGAC));
                colorbarXRangeGAC = 128.5 + 127.5 * colorbarXPercentGAC/100;
                colorbarXColorsGAC = round(256*(interp1(1:256,colorMap,colorbarXRangeGAC)));
                for iYPos = 1: numel(colorbarYPosGAC)
                    showImageGAC(colorbarYPosGAC(iYPos),colorbarXPosGAC,:) = colorbarXColorsGAC;
                end
                % percent scale
                colorbarXPercentLabelsGAC = 0:10:100;
                colorbarXPercentLabelsPosGAC = colorbarXPosGAC(arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXPercentLabelsGAC));
                showImageGAC((colorbarYPosGAC(end)+1):(colorbarYPosGAC(end)+5),colorbarXPercentLabelsPosGAC,:) = 255;

                % subjects scale
                stepsizeGAC = ceil(numel(subjectsAboveCriterion)/10);
                colorbarXSubjectLabelsGAC = 1:stepsizeGAC:numel(subjectsAboveCriterion);
                colorbarXSubjectLabelsGACfine = 1:numel(subjectsAboveCriterion);
                colorbarXSubjectPercentLabelsGAC = 100*colorbarXSubjectLabelsGAC/numel(subjectsAboveCriterion);
                colorbarXSubjectPercentLabelsGACfine = 100*colorbarXSubjectLabelsGACfine/numel(subjectsAboveCriterion);
                temp = arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXSubjectPercentLabelsGAC, 'UniformOutput', false);
                tempfine = arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXSubjectPercentLabelsGACfine, 'UniformOutput', false);
                colorbarXSubjectPercentLabelsPosGAC = colorbarXPosGAC(cellfun(@(x) x(1), temp));
                colorbarXSubjectPercentLabelsPosGACfine = colorbarXPosGAC(cellfun(@(x) x(1), tempfine));

                showImageGAC((colorbarYPosGAC(1)-5):(colorbarYPosGAC(1)-1),colorbarXSubjectPercentLabelsPosGAC,:) = 255;
                showImageGAC(colorbarYPosGAC,(colorbarXPosGAC(1)):(colorbarXSubjectPercentLabelsPosGAC(1)-1),:) = 0;
                showImageGAC(colorbarYPosGAC,(colorbarXSubjectPercentLabelsPosGACfine(max(unique(gac)))+1):(colorbarXPosGAC(end)),:) = 0;
            end

            % difference map
            showImageGBCminusGAC = uint8(zeros(size(mni,1), size(mni,2), 3));
            showImageGBCminusGAC(:,:,1) = squeeze(redGBCminusGAC(:,:,useSlice));
            showImageGBCminusGAC(:,:,2) = squeeze(greenGBCminusGAC(:,:,useSlice));
            showImageGBCminusGAC(:,:,3) = squeeze(blueGBCminusGAC(:,:,useSlice));
            showImageGBCminusGAC(size(mni,1)+80, size(mni,2),3) = 0;
            % colorbar
            colorbarYPosGBCminusGAC = size(mni,1)+10:size(mni,1)+20;
            colorbarXPosGBCminusGAC = 20:size(mni,2)-20;
            showImageGBCminusGAC((colorbarYPosGBCminusGAC(1)-1):(colorbarYPosGBCminusGAC(end)+1),(colorbarXPosGBCminusGAC(1)-1):(colorbarXPosGBCminusGAC(end)+1),:) = 255;
            % map to [-1:1]  round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),1));
            colorbarXPercentGBCminusGAC = 200*(colorbarXPosGBCminusGAC-min(colorbarXPosGBCminusGAC))/(max(colorbarXPosGBCminusGAC)-min(colorbarXPosGBCminusGAC))-100;
            colorbarXRangeGBCminusGAC = 128.5 + 127.5 * colorbarXPercentGBCminusGAC/100;
            colorbarXColorsGBCminusGAC = round(256*(interp1(1:256,colorMap,colorbarXRangeGBCminusGAC)));
            for iYPos = 1: numel(colorbarYPosGBCminusGAC)
                showImageGBCminusGAC(colorbarYPosGBCminusGAC(iYPos),colorbarXPosGBCminusGAC,:) = colorbarXColorsGBCminusGAC;
            end
            % percent scale
            colorbarXPercentLabelsGBCminusGAC = -100:10:100;
            colorbarXPercentLabelsPosGBCminusGAC = colorbarXPosGBCminusGAC(arrayfun(@(x) find(min(abs(colorbarXPercentGBCminusGAC-x)) == abs(colorbarXPercentGBCminusGAC-x)), colorbarXPercentLabelsGBCminusGAC));
            showImageGBCminusGAC((colorbarYPosGBCminusGAC(end)+1):(colorbarYPosGBCminusGAC(end)+5),colorbarXPercentLabelsPosGBCminusGAC,:) = 255;

            importantPoints = 100 * [min(uniquegbcminusgac), -params.cutoffDifferenceImage, params.cutoffDifferenceImage, max(uniquegbcminusgac)];
            temp = arrayfun(@(x) find(min(abs(colorbarXPercentGBCminusGAC-x)) == abs(colorbarXPercentGBCminusGAC-x)), importantPoints, 'UniformOutput', false);
            importantPointsX = colorbarXPosGBCminusGAC(cellfun(@(x) x(1), temp));


            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(colorbarXPosGBCminusGAC(1)):(importantPointsX(1)-1),:) = 0;
            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(importantPointsX(2)+1):(importantPointsX(3)-1),:) = 0;
            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(importantPointsX(end)+1):(colorbarXPosGBCminusGAC(end)),:) = 0;





            showImageGBCminusGACDetail = showImageGBCminusGAC;
            % borders of ROIs

            for iROI = 1: numel(params.showRegions)
                ROISlice = logical(squeeze(ROI{iROI,1}(:,:,useSlice)));
                L = bwboundaries(ROISlice);
                for iL = 1:numel(L)
                    for iPix = 1: size(L{iL},1)
                        showImageGBC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                        showImageGAC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                        showImageGBCminusGAC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;

                    end
                end
            end


            % borders of ROIs for Details plot
            if ~params.showOnlyRegionsDetails
                showImageGBCminusGACDetail = showImageGBCminusGAC;
            end
            for iROI = 1: numel(params.showRegionsDetail)
                ROISliceDetail = logical(squeeze(ROIDetail{iROI,1}(:,:,useSlice)));
                L = bwboundaries(ROISliceDetail);
                RoiDetailUsed(iROI) = ~isempty(L);
                for iL = 1:numel(L)
                    for iPix = 1: size(L{iL},1)
                        if params.RegionDetailsMonochrome
                            showImageGBCminusGACDetail(L{iL}(iPix,1),L{iL}(iPix,2),:)=[1 1 1]*params.ROIColors;
                        else
                            showImageGBCminusGACDetail(L{iL}(iPix,1),L{iL}(iPix,2),:)=textcol{rem(iROI-1,nROIsDetail)+1}*255;
                        end
                    end
                end
            end
            % details
            % difference map
            showDetailsGBCminusGAC = uint8(zeros(size(mni,1), size(mni,2), 3));
            xstart = size(mni,1)/4;
            xend = xstart + size(mni,1)/2 -1;
            ystart = size(mni,2)/4;
            yend = ystart + size(mni,2)/2 -1;

            size(showDetailsGBCminusGAC)
            showDetailsGBCminusGAC(1:2:end,1:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(2:2:end,1:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(1:2:end,2:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(2:2:end,2:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);

            im1 = permute(showImageGBC, [2,1,3]);
            im2 = permute(showImageGAC, [2,1,3]);
            im3 = permute(showImageGBCminusGAC, [2,1,3]);
            im4 = permute(showDetailsGBCminusGAC, [2,1,3]);

            %% axial figure
            fig = figure;
            fig.WindowState = 'maximized';
            im = [im1(end:-1:1,:,:), im2(end:-1:1,:,:), im3(end:-1:1,:,:), im4(end:-1:1,:,:)];
            imshow(im);
            if params.doFlipping
                text(round(size(im,2)*0.75)+70, 20,'IA', 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment', 'left')
            else
                text(round(size(im,2)*0.75)+70, 20,'RA', 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment', 'left')
            end
            text(size(im,2)-20, 20,sprintf('z = %i', round(MNI_Z)), 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment', 'right')
            % GBC
            for iLabel = 1:2:numel(colorbarXPercentLabelsGBC)
                text(colorbarYPosGBC(end)+10, size(im,1) - colorbarXPercentLabelsPosGBC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
            end

            for iLabel = 1: numel(colorbarXSubjectPercentLabelsGBC)
                if colorbarXSubjectLabelsGBC(iLabel) <= max(unique(gbc))
                    text(colorbarYPosGBC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGBC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [1,1,1], 'HorizontalAlignment','right')
                else
                    text(colorbarYPosGBC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGBC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment','right')
                end

            end
            % GAC
            if numel(subjectsAboveCriterion)>0
                for iLabel = 1:2:numel(colorbarXPercentLabelsGAC)
                    text(size(mni,1) + 50 + colorbarYPosGAC(end)+10,  size(im,1) - colorbarXPercentLabelsPosGAC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
                end
                tempMax = max(unique(gac));
                for iLabel = 1: numel(colorbarXSubjectPercentLabelsGAC)
                    if colorbarXSubjectLabelsGAC(iLabel) <= tempMax
                        text(size(mni,1) + 50 + colorbarYPosGAC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGAC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1], 'HorizontalAlignment','right')
                    else
                        text(size(mni,1) + 50 + colorbarYPosGAC(1)-10,  size(im,1) - colorbarXSubjectPercentLabelsPosGAC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment','right')
                    end

                end
            end
            % GBC-GAC
            for iLabel = 1:2:numel(colorbarXPercentLabelsGBCminusGAC)
                text(2* size(mni,1) + 100 + colorbarYPosGBCminusGAC(end)+10,  size(im,1) - colorbarXPercentLabelsPosGBCminusGAC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGBCminusGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
            end
            %ROI
            if ~params.RegionDetailsMonochrome
                lineText = 1;
                maxTextLength = max(cellfun(@numel,textstr));
                RoiDetailUsed = RoiDetailUsed(1:numel(RoiDetailUsed)/2) | RoiDetailUsed((numel(RoiDetailUsed)/2)+1:end);
                for iROI = find(RoiDetailUsed)
                    while numel(textstr{iROI}) < maxTextLength
                        textstr{iROI} = [textstr{iROI}, ' '];
                    end
                    text(round(size(im,2)*0.75)+70, 10 + 60 + lineText * 24, textstr{iROI}, 'Color', [1,1,1], 'FontName', 'FixedWidth', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none', 'BackgroundColor', textcol{iROI})
                    lineText = lineText + 1;
                end
            end

            % generate saveName and save figure
            if iDeficit == 5
                saveName = [figFolder, filesep, sprintf('LesionOverlap_z_%0.2f_slice_%03.0f_combined', MNI_Z, useSlice)];
            else
                saveName = [figFolder, filesep, sprintf('%s_threshold_%0.1f_z_%0.2f_slice_%03.0f_combined',params.criterionVariable, params.criterionThreshold, MNI_Z, useSlice)];
            end

            if params.doThalamusMasking
                saveName = [saveName, '_masked'];
            end
            if params.doFlipping
                saveName = [saveName, '_flipped'];
            end
            if params.replaceColorMap
                saveName = [saveName, sprintf('_ColorSaturationAt%0.2f', params.colorSaturationAt)];
            end
            if params.RegionDetailsMonochrome
                saveName = [saveName, '_ROIDetailsMonoChrome']
            end
            %saveName = [saveName, '_', sprintf('_%s', string(datetime('now','format','yyyy_MM_dd_HH_mm')))];
            axis equal
            axis tight
            axis off
            truesize
            savefig(fig, [saveName, '.fig']);

        end

   
        %% plotSlices Coronal
        for useSlice =  params.plotSlicesY
            MNI_Y = mniHeader.Transform.T(4,2) + (useSlice-0.5)* mniHeader.Transform.T(2,2);
            showImageGBC = uint8(zeros(size(mni,1), size(mni,3), 3));
            showImageGBC(:,:,1) = squeeze(redGBC(:,useSlice,:));
            showImageGBC(:,:,2) = squeeze(greenGBC(:,useSlice,:));
            showImageGBC(:,:,3) = squeeze(blueGBC(:,useSlice,:));
            showImageGBC(size(mni,1)+50, size(mni,3),3) = 0;
            % colorbar
            colorbarYPosGBC = size(mni,1)+10:size(mni,1)+20;
            colorbarXPosGBC = size(mni,3)/2:size(mni,3)-20;
            showImageGBC((colorbarYPosGBC(1)-1):(colorbarYPosGBC(end)+1),(colorbarXPosGBC(1)-1):(colorbarXPosGBC(end)+1),:) = 255;

            % map to [0:1]
            colorbarXPercentGBC = 100*(colorbarXPosGBC-min(colorbarXPosGBC))/(max(colorbarXPosGBC)-min(colorbarXPosGBC));
            colorbarXRangeGBC = 128.5 + 127.5 * colorbarXPercentGBC/100;
            colorbarXColorsGBC = round(256*(interp1(1:256,colorMap,colorbarXRangeGBC)));
            for iYPos = 1: numel(colorbarYPosGBC)
                showImageGBC(colorbarYPosGBC(iYPos),colorbarXPosGBC,:) = colorbarXColorsGBC;
            end
            % percent scale
            colorbarXPercentLabelsGBC = 0:10:100;
            colorbarXPercentLabelsPosGBC = colorbarXPosGBC(arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXPercentLabelsGBC));
            showImageGBC((colorbarYPosGBC(end)+1):(colorbarYPosGBC(end)+5),colorbarXPercentLabelsPosGBC,:) = 255;
            % subjects scale
            stepsizeGBC = ceil(numel(subjectsBelowCriterion)/10);
            colorbarXSubjectLabelsGBC = 1:stepsizeGBC:numel(subjectsBelowCriterion);
            colorbarXSubjectLabelsGBCfine = 1:numel(subjectsBelowCriterion)
            colorbarXSubjectPercentLabelsGBC = 100*colorbarXSubjectLabelsGBC/numel(subjectsBelowCriterion);
            colorbarXSubjectPercentLabelsGBCfine = 100*colorbarXSubjectLabelsGBCfine/numel(subjectsBelowCriterion);
            temp = arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXSubjectPercentLabelsGBC, 'UniformOutput', false);
            tempfine = arrayfun(@(x) find(min(abs(colorbarXPercentGBC-x)) == abs(colorbarXPercentGBC-x)), colorbarXSubjectPercentLabelsGBCfine, 'UniformOutput', false);
            colorbarXSubjectPercentLabelsPosGBC = colorbarXPosGBC(cellfun(@(x) x(1), temp));
            colorbarXSubjectPercentLabelsPosGBCfine = colorbarXPosGBC(cellfun(@(x) x(1), tempfine));

            showImageGBC((colorbarYPosGBC(1)-5):(colorbarYPosGBC(1)-1),colorbarXSubjectPercentLabelsPosGBC,:) = 255;
            showImageGBC(colorbarYPosGBC,(colorbarXPosGBC(1)):(colorbarXSubjectPercentLabelsPosGBC(1)-1),:) = 0;
            showImageGBC(colorbarYPosGBC,(colorbarXSubjectPercentLabelsPosGBCfine(max(unique(gbc)))+1):(colorbarXPosGBC(end)),:) = 0;



            showImageGAC = uint8(zeros(size(mni,1), size(mni,3), 3));
            showImageGAC(:,:,1) = squeeze(redGAC(:,useSlice,:));
            showImageGAC(:,:,2) = squeeze(greenGAC(:,useSlice,:));
            showImageGAC(:,:,3) = squeeze(blueGAC(:,useSlice,:));
            showImageGAC(size(mni,1)+50, size(mni,3),3) = 0;
            if numel(subjectsAboveCriterion) > 0

                % colorbar
                colorbarYPosGAC = size(mni,1)+10:size(mni,1)+20;
                colorbarXPosGAC = size(mni,3)/2:size(mni,3)-20;
                showImageGAC((colorbarYPosGAC(1)-1):(colorbarYPosGAC(end)+1),(colorbarXPosGAC(1)-1):(colorbarXPosGAC(end)+1),:) = 255;

                % map to [0:1]
                colorbarXPercentGAC = 100*(colorbarXPosGAC-min(colorbarXPosGAC))/(max(colorbarXPosGAC)-min(colorbarXPosGAC));
                colorbarXRangeGAC = 128.5 + 127.5 * colorbarXPercentGAC/100;
                colorbarXColorsGAC = round(256*(interp1(1:256,colorMap,colorbarXRangeGAC)));
                for iYPos = 1: numel(colorbarYPosGAC)
                    showImageGAC(colorbarYPosGAC(iYPos),colorbarXPosGAC,:) = colorbarXColorsGAC;
                end
                % percent scale
                colorbarXPercentLabelsGAC = 0:10:100;
                colorbarXPercentLabelsPosGAC = colorbarXPosGAC(arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXPercentLabelsGAC));
                showImageGAC((colorbarYPosGAC(end)+1):(colorbarYPosGAC(end)+5),colorbarXPercentLabelsPosGAC,:) = 255;
                % subjects scale
                stepsizeGAC = ceil(numel(subjectsAboveCriterion)/10);
                colorbarXSubjectLabelsGAC = 1:stepsizeGAC:numel(subjectsAboveCriterion);
                colorbarXSubjectLabelsGACfine = 1:numel(subjectsAboveCriterion);
                colorbarXSubjectPercentLabelsGAC = 100*colorbarXSubjectLabelsGAC/numel(subjectsAboveCriterion);
                colorbarXSubjectPercentLabelsGACfine = 100*colorbarXSubjectLabelsGACfine/numel(subjectsAboveCriterion);
                temp = arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXSubjectPercentLabelsGAC, 'UniformOutput', false);
                tempfine = arrayfun(@(x) find(min(abs(colorbarXPercentGAC-x)) == abs(colorbarXPercentGAC-x)), colorbarXSubjectPercentLabelsGACfine, 'UniformOutput', false);
                colorbarXSubjectPercentLabelsPosGAC = colorbarXPosGAC(cellfun(@(x) x(1), temp));
                colorbarXSubjectPercentLabelsPosGACfine = colorbarXPosGAC(cellfun(@(x) x(1), tempfine));

                showImageGAC((colorbarYPosGAC(1)-5):(colorbarYPosGAC(1)-1),colorbarXSubjectPercentLabelsPosGAC,:) = 255;
                showImageGAC(colorbarYPosGAC,(colorbarXPosGAC(1)):(colorbarXSubjectPercentLabelsPosGAC(1)-1),:) = 0;
                showImageGAC(colorbarYPosGAC,(colorbarXSubjectPercentLabelsPosGACfine(max(unique(gac)))+1):(colorbarXPosGAC(end)),:) = 0;

            end

            % difference map
            showImageGBCminusGAC = uint8(zeros(size(mni,1), size(mni,3), 3));
            showImageGBCminusGAC(:,:,1) = squeeze(redGBCminusGAC(:,useSlice,:));
            showImageGBCminusGAC(:,:,2) = squeeze(greenGBCminusGAC(:,useSlice,:));
            showImageGBCminusGAC(:,:,3) = squeeze(blueGBCminusGAC(:,useSlice,:));
            showImageGBCminusGAC(size(mni,1)+80, size(mni,3),3) = 0;
            % colorbar
            colorbarYPosGBCminusGAC = size(mni,1)+10:size(mni,1)+20;
            colorbarXPosGBCminusGAC = 20:size(mni,3)-20;
            showImageGBCminusGAC((colorbarYPosGBCminusGAC(1)-1):(colorbarYPosGBCminusGAC(end)+1),(colorbarXPosGBCminusGAC(1)-1):(colorbarXPosGBCminusGAC(end)+1),:) = 255;
            % map to [-1:1]  round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),1));
            colorbarXPercentGBCminusGAC = 200*(colorbarXPosGBCminusGAC-min(colorbarXPosGBCminusGAC))/(max(colorbarXPosGBCminusGAC)-min(colorbarXPosGBCminusGAC))-100;
            colorbarXRangeGBCminusGAC = 128.5 + 127.5 * colorbarXPercentGBCminusGAC/100;
            colorbarXColorsGBCminusGAC = round(256*(interp1(1:256,colorMap,colorbarXRangeGBCminusGAC)));
            for iYPos = 1: numel(colorbarYPosGBCminusGAC)
                showImageGBCminusGAC(colorbarYPosGBCminusGAC(iYPos),colorbarXPosGBCminusGAC,:) = colorbarXColorsGBCminusGAC;
            end
            % percent scale
            colorbarXPercentLabelsGBCminusGAC = -100:10:100;
            colorbarXPercentLabelsPosGBCminusGAC = colorbarXPosGBCminusGAC(arrayfun(@(x) find(min(abs(colorbarXPercentGBCminusGAC-x)) == abs(colorbarXPercentGBCminusGAC-x)), colorbarXPercentLabelsGBCminusGAC));
            showImageGBCminusGAC((colorbarYPosGBCminusGAC(end)+1):(colorbarYPosGBCminusGAC(end)+5),colorbarXPercentLabelsPosGBCminusGAC,:) = 255;

            importantPoints = 100 * [min(uniquegbcminusgac), -params.cutoffDifferenceImage, params.cutoffDifferenceImage, max(uniquegbcminusgac)];
            temp = arrayfun(@(x) find(min(abs(colorbarXPercentGBCminusGAC-x)) == abs(colorbarXPercentGBCminusGAC-x)), importantPoints, 'UniformOutput', false);
            importantPointsX = colorbarXPosGBCminusGAC(cellfun(@(x) x(1), temp));


            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(colorbarXPosGBCminusGAC(1)):(importantPointsX(1)-1),:) = 0;
            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(importantPointsX(2)+1):(importantPointsX(3)-1),:) = 0;
            showImageGBCminusGAC(colorbarYPosGBCminusGAC,(importantPointsX(end)+1):(colorbarXPosGBCminusGAC(end)),:) = 0;
            showImageGBCminusGACDetail = showImageGBCminusGAC;
            % borders of ROIs
            for iROI = 1: numel(params.showRegions)
                ROISlice = logical(squeeze(ROI{iROI,1}(:,useSlice,:)));
                L = bwboundaries(ROISlice);
                for iL = 1:numel(L)
                    for iPix = 1: size(L{iL},1)
                        showImageGBC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                        showImageGAC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                        showImageGBCminusGAC(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;

                    end
                end
            end
            % borders of ROIs for Details plot
            if ~params.showOnlyRegionsDetails
                showImageGBCminusGACDetail = showImageGBCminusGAC;
            end

            for iROI = 1: numel(params.showRegionsDetail)
                ROISliceDetail = logical(squeeze(ROIDetail{iROI,1}(:,useSlice,:)));
                L = bwboundaries(ROISliceDetail);
                RoiDetailUsed(iROI) = ~isempty(L);
                for iL = 1:numel(L)

                    for iPix = 1: size(L{iL},1)
                        if params.RegionDetailsMonochrome
                            showImageGBCminusGACDetail(L{iL}(iPix,1),L{iL}(iPix,2),:)=[1 1 1]*params.ROIColors;
                        else
                            showImageGBCminusGACDetail(L{iL}(iPix,1),L{iL}(iPix,2),:)=textcol{rem(iROI-1,nROIsDetail)+1}*255;
                        end
                    end
                end
            end

            % details
            % difference map
            showDetailsGBCminusGAC = uint8(zeros(size(mni,1), size(mni,3), 3));
            xstart = size(mni,1)/4;
            xend = xstart + size(mni,1)/2 -1;
            ystart = size(mni,3)/4;
            yend = ystart + size(mni,3)/2 -1;

            size(showDetailsGBCminusGAC)
            showDetailsGBCminusGAC(1:2:end,1:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(2:2:end,1:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(1:2:end,2:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);
            showDetailsGBCminusGAC(2:2:end,2:2:end,:) = showImageGBCminusGACDetail(xstart:xend,ystart:yend,:);

            im1 = permute(showImageGBC, [2,1,3]);
            im2 = permute(showImageGAC, [2,1,3]);
            im3 = permute(showImageGBCminusGAC, [2,1,3]);
            im4 = permute(showDetailsGBCminusGAC, [2,1,3]);

            %% slice figure
            fig = figure;
            fig.WindowState = 'maximized';
            im = [im1(end:-1:1,:,:), im2(end:-1:1,:,:), im3(end:-1:1,:,:), im4(end:-1:1,:,:)];
            imshow(im);
            if params.doFlipping
                text(round(size(im,2)*0.75)+70, 20,'IS', 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment','left')
            else
                text(round(size(im,2)*0.75)+70, 20,'RS', 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment','left')
            end
            text(size(im,2)-20, 20,sprintf('y = %i', round(MNI_Y)), 'FontSize', 16, 'Color', [1,1,1], 'HorizontalAlignment','right')
            % GBC
            for iLabel = 1:2:numel(colorbarXPercentLabelsGBC)
                text(colorbarYPosGBC(end)+10, size(im,1) - colorbarXPercentLabelsPosGBC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
            end
            for iLabel = 1: numel(colorbarXSubjectPercentLabelsGBC)
                if colorbarXSubjectLabelsGBC(iLabel) <= max(unique(gbc))
                    text(colorbarYPosGBC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGBC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [1,1,1], 'HorizontalAlignment','right')
                else
                    text(colorbarYPosGBC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGBC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGBC(iLabel)), 'FontSize', 12, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment','right')
                end

            end
            % GAC
            if numel(subjectsAboveCriterion)>0
                for iLabel = 1:2:numel(colorbarXPercentLabelsGAC)
                    text(size(mni,1) + 50 + colorbarYPosGAC(end)+10,  size(im,1) - colorbarXPercentLabelsPosGAC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
                end
                tempMax = max(unique(gac));
                for iLabel = 1: numel(colorbarXSubjectPercentLabelsGAC)
                    if colorbarXSubjectLabelsGAC(iLabel) <= tempMax
                        text(size(mni,1) + 50 + colorbarYPosGAC(1)-10, size(im,1) - colorbarXSubjectPercentLabelsPosGAC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1], 'HorizontalAlignment','right')
                    else
                        text(size(mni,1) + 50 + colorbarYPosGAC(1)-10,  size(im,1) - colorbarXSubjectPercentLabelsPosGAC(iLabel),  sprintf('%i', colorbarXSubjectLabelsGAC(iLabel)), 'FontSize', 12, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment','right')
                    end

                end
            end
            % GBC-GAC
            for iLabel = 1:2:numel(colorbarXPercentLabelsGBCminusGAC)
                text(2* size(mni,1) + 100 + colorbarYPosGBCminusGAC(end)+10,  size(im,1) - colorbarXPercentLabelsPosGBCminusGAC(iLabel),  sprintf('%i %%', colorbarXPercentLabelsGBCminusGAC(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
            end
            if ~params.RegionDetailsMonochrome
                %ROI
                lineText = 1;
                maxTextLength = max(cellfun(@numel,textstr));
                RoiDetailUsed = RoiDetailUsed(1:numel(RoiDetailUsed)/2) | RoiDetailUsed((numel(RoiDetailUsed)/2)+1:end);
                for iROI = find(RoiDetailUsed)
                    while numel(textstr{iROI}) < maxTextLength
                        textstr{iROI} = [textstr{iROI}, ' '];
                    end
                    text(round(size(im,2)*0.75)+70, 10 + 60 + lineText * 24, textstr{iROI}, 'Color', [1,1,1], 'FontName', 'FixedWidth', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none', 'BackgroundColor', textcol{iROI})
                    lineText = lineText + 1;
                end

            end

            % generate saveName and save figure
            if iDeficit == 5
                saveName = [figFolder, filesep, sprintf('LesionOverlap_y_%0.2f_slice_%03.0f_combined', MNI_Y, useSlice)];
            else
                saveName = [figFolder, filesep, sprintf('%s_threshold_%0.1f_y_%0.2f_slice_%03.0f_combined',params.criterionVariable, params.criterionThreshold, MNI_Y, useSlice)];
            end
            if params.doThalamusMasking
                saveName = [saveName, '_masked'];
            end
            if params.doFlipping
                saveName = [saveName, '_flipped'];
            end
            if params.replaceColorMap
                saveName = [saveName, sprintf('_ColorSaturationAt%0.2f', params.colorSaturationAt)];
            end
            if params.RegionDetailsMonochrome
                saveName = [saveName, '_ROIDetailsMonoChrome']
            end
            %saveName = [saveName, '_', sprintf('_%s', string(datetime('now','format','yyyy_MM_dd_HH_mm')))];
          
            axis equal
            axis tight
            axis off
            truesize
            savefig(fig, [saveName, '.fig']);
        end


        %% combined plots
        for iVisual = 4 %1:4
            switch iVisual
                case 1
                    params.RegionDetailsMonochrome = false;
                    params.showOnlyRegionsDetails = false;

                case 2
                    params.RegionDetailsMonochrome = true;
                    params.showOnlyRegionsDetails = false;
                case 3
                    params.RegionDetailsMonochrome = false;
                    params.showOnlyRegionsDetails = true;
                case 4
                    params.RegionDetailsMonochrome = true;
                    params.showOnlyRegionsDetails = true;

            end
            layout.data = cell(2,6);
            layout.orientation = [2*ones(1,6); 3*ones(1,6)];
            layout.slice = [224, 220, 216, 212, 208, 212;...
                150, 154, 158, 162, 166, 162];
            layout.doEnlarge = [false(2,5), true(2,1)];

            layout.slice = [224, 220, 216, 212, 208, 204;...
                146, 150, 154, 158, 162, 166];
            layout.doEnlarge = [ true(2,6)];
            layout.slice = [226, 222, 218, 214, 210, 206;...
                150, 154, 158, 162, 166, 170];
            layout.doEnlarge = [ true(2,6)];
 
            % changed to 5 images to fit the screen for images in original
            % size
             layout.data = cell(2,5);
             layout.slice = [226, 220, 214, 208, 202;...
                 148, 154, 160, 166, 172];
             layout.doEnlarge = [ true(2,5)];

            % recalculate difference image with params.cutoffDifferenceImage
            redGBCminusGAC = mni;
            greenGBCminusGAC = mni;
            blueGBCminusGAC =  mni;
            uniquegbcminusgac = unique(gbcminusgac);
            for i = 1:numel(uniquegbcminusgac)
                %uniquegbcminusgac(i)
                if abs(uniquegbcminusgac(i)) >= params.cutoffDifferenceImage
                    redGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),1));
                    greenGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),2));
                    blueGBCminusGAC(gbcminusgac==uniquegbcminusgac(i)) = round(256 * colorMap(round(128.5 + 127.5*uniquegbcminusgac(i)),3));
                end
            end


            for r = 1:size(layout.data,1)
                for c = 1: size(layout.data,2)
                    switch layout.orientation(r,c)
                        case 1
                            % todo: prepare for axial plots
                            % PROB: different sizes of images
                            layout.MNI(r,c) = round(mniHeader.Transform.T(4,1) + (layout.slice(r,c)-0.5)* mniHeader.Transform.T(1,1));
                            layout.sliceLabel{r,c} = sprintf('x = %i', layout.MNI(r,c));

                        case 2
                            layout.MNI(r,c) = round(mniHeader.Transform.T(4,2) + (layout.slice(r,c)-0.5)* mniHeader.Transform.T(2,2));
                            layout.sliceLabel{r,c} = sprintf('y = %i', layout.MNI(r,c));
                            if params.doFlipping
                                layout.orientationLabel{r,c} = 'IS';
                            else
                                layout.orientationLabel{r,c} = 'RS';
                            end
                            layout.data{r,c} = uint8(zeros(size(mni,1), size(mni,3), 3));
                            layout.data{r,c}(:,:,1) = squeeze(redGBCminusGAC(:,layout.slice(r,c),:));
                            layout.data{r,c}(:,:,2) = squeeze(greenGBCminusGAC(:,layout.slice(r,c),:));
                            layout.data{r,c}(:,:,3) = squeeze(blueGBCminusGAC(:,layout.slice(r,c),:));
                            % borders of ROIs
                            if ~layout.doEnlarge(r,c) || ~params.showOnlyRegionsDetails
                                for iROI = 1: numel(params.showRegions)
                                    ROISlice = logical(squeeze(ROI{iROI,1}(:,layout.slice(r,c),:)));
                                    L = bwboundaries(ROISlice);
                                    for iL = 1:numel(L)
                                        for iPix = 1: size(L{iL},1)
                                            layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                                        end
                                    end
                                end
                            end
                            % borders of ROIs for Details plot
                            if layout.doEnlarge(r,c)
                                for iROI = 1: numel(params.showRegionsDetail)
                                    ROISliceDetail = logical(squeeze(ROIDetail{iROI,1}(:,layout.slice(r,c),:)));
                                    L = bwboundaries(ROISliceDetail);
                                    layout.RoiDetailUsed{r,c}(iROI) = ~isempty(L);
                                    for iL = 1:numel(L)
                                        for iPix = 1: size(L{iL},1)
                                            if params.RegionDetailsMonochrome
                                                layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                                            else
                                                layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=textcol{rem(iROI-1,nROIsDetail)+1}*255;
                                            end
                                        end
                                    end
                                end
                            end

                        case 3
                            layout.MNI(r,c) = round(mniHeader.Transform.T(4,3) + (layout.slice(r,c)-0.5)* mniHeader.Transform.T(3,3));
                            layout.sliceLabel{r,c} = sprintf('z = %i', layout.MNI(r,c));
                            if params.doFlipping
                                layout.orientationLabel{r,c} = 'IA';
                            else
                                layout.orientationLabel{r,c} = 'RA';
                            end
                            layout.data{r,c} = uint8(zeros(size(mni,1), size(mni,2), 3));
                            layout.data{r,c}(:,:,1) = squeeze(redGBCminusGAC(:,:,layout.slice(r,c)));
                            layout.data{r,c}(:,:,2) = squeeze(greenGBCminusGAC(:,:,layout.slice(r,c)));
                            layout.data{r,c}(:,:,3) = squeeze(blueGBCminusGAC(:,:,layout.slice(r,c)));
                            % borders of ROIs
                            if ~layout.doEnlarge(r,c) || ~params.showOnlyRegionsDetails
                                for iROI = 1: numel(params.showRegions)
                                    ROISlice = logical(squeeze(ROI{iROI,1}(:,:,layout.slice(r,c))));
                                    L = bwboundaries(ROISlice);
                                    for iL = 1:numel(L)
                                        for iPix = 1: size(L{iL},1)
                                            layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                                        end
                                    end
                                end
                            end
                            % borders of ROIs for Details plot
                            if layout.doEnlarge(r,c)
                                for iROI = 1: numel(params.showRegionsDetail)
                                    ROISliceDetail = logical(squeeze(ROIDetail{iROI,1}(:,:,layout.slice(r,c))));
                                    L = bwboundaries(ROISliceDetail);
                                    layout.RoiDetailUsed{r,c}(iROI) = ~isempty(L);
                                    for iL = 1:numel(L)
                                        for iPix = 1: size(L{iL},1)
                                            if params.RegionDetailsMonochrome
                                                layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=params.ROIColors;
                                            else
                                                layout.data{r,c}(L{iL}(iPix,1),L{iL}(iPix,2),:)=textcol{rem(iROI-1,nROIsDetail)+1}*255;
                                            end
                                        end
                                    end
                                end
                            end
                        otherwise
                            error('unknown orientation')
                    end

                    if layout.doEnlarge(r,c)
                        imagesize = size(layout.data{r,c});
                        imageEnlarged = layout.data{r,c}(imagesize(1)/4 : (imagesize(1)*3/4) - 1, imagesize(2)/4 : (imagesize(2)*3/4) - 1, :);
                        layout.data{r,c}(1:2:end,1:2:end,:) = imageEnlarged;
                        layout.data{r,c}(2:2:end,1:2:end,:) = imageEnlarged;
                        layout.data{r,c}(1:2:end,2:2:end,:) = imageEnlarged;
                        layout.data{r,c}(2:2:end,2:2:end,:) = imageEnlarged;
                    end

                    layout.data{r,c} = permute(layout.data{r,c},[2,1,3]);
                    layout.data{r,c} = layout.data{r,c}(end:-1:1,:,:);

                end
            end

            canvas = [];
            for r = 1:size(layout.data,1)
                rowcanvas = [];
                for c = 1: size(layout.data,2)
                    layout.position{r,c} = [size(rowcanvas,2) + 1, size(canvas,1) + 1, size(layout.data{r,c}, 2), size(layout.data{r,c}, 1)];
                    rowcanvas = [rowcanvas, layout.data{r,c}];

                end
                canvas =[canvas;rowcanvas];
            end

            % add colorbar
            % canvas x will become rows, y will become colums

            % black area on the right
            canvas(size(canvas,1), size(canvas,2) + 180, :) = 0;
            % horizontal extent of the color bar (= second coordinate (y) is columns)
            colorbarYPosLesion = size(canvas,2)-120:size(canvas,2)-80;
            % vertical extent is first coordinate (x) is rows.
            colorbarXPosLesion = 100:size(canvas,1)-100;
            % precolor area in white (draws white borders)
            canvas((colorbarXPosLesion(1)-1):(colorbarXPosLesion(end)+1),(colorbarYPosLesion(1)-1):(colorbarYPosLesion(end)+1),:) = 255;
            if iDeficit == 5
                % map color bar to [0:1] ([0%:100%])
                %100*(colorbarXPosGBCminusGAC-min(colorbarXPosGBCminusGAC))/(max(colorbarXPosGBCminusGAC)-min(colorbarXPosGBCminusGAC));
                colorbarXPercentLesion = 100*(colorbarXPosLesion-min(colorbarXPosLesion))/(max(colorbarXPosLesion)-min(colorbarXPosLesion));
                colorbarXRangeLesion = 128.5 + 127.5 * colorbarXPercentLesion/100;
                colorbarXColorsLesion = round(256*(interp1(1:256,colorMap,colorbarXRangeLesion)));
                % paint the colorbar
                for iYPos = 1: numel(colorbarYPosLesion)
                    canvas(size(canvas,1)-colorbarXPosLesion,colorbarYPosLesion(iYPos),:) = colorbarXColorsLesion;
                end
                % percent scale
                colorbarXPercentLabelsLesion = 0:10:100;
                colorbarXPercentLabelsPosLesion = colorbarXPosLesion(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), colorbarXPercentLabelsLesion));
                % add strokes to scale
                canvas( colorbarXPercentLabelsPosLesion, (colorbarYPosLesion(end)+1):(colorbarYPosLesion(end)+5), :) = 255;
    
                % subjects scale (1:nSubj)
                colorbarXSubjectLabelsLesion = 0:numel(lesionFilenames);
                colorbarXSubjectPercentLabelsLesion = 100*colorbarXSubjectLabelsLesion/numel(lesionFilenames);
                temp = arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), colorbarXSubjectPercentLabelsLesion, 'UniformOutput', false);
                colorbarXSubjectPercentLabelsPosLesion = colorbarXPosLesion(cellfun(@(x) x(1), temp));
                % add strokes to scale
                canvas(colorbarXSubjectPercentLabelsPosLesion,(colorbarYPosLesion(1)-5):(colorbarYPosLesion(1)-1),:) = 255;
                 
    

                % remove color from scale of less than abs(params.cutoffDifferenceImage)
                lowerCutoff = -params.cutoffDifferenceImage * 100;
                upperCutoff = params.cutoffDifferenceImage * 100;
                mask = [lowerCutoff, upperCutoff];
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), mask, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(temp(1))):(colorbarXPosLesion(temp(2))-1)),colorbarYPosLesion,:) = 0;
                % remove color from scale above max percentage
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), max(uniquegbcminusgac)*100, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(temp):colorbarXPosLesion(end)-1)),colorbarYPosLesion,:) = 0;
                % remove color from scale below min percentage
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), min(uniquegbcminusgac)*100, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(1):colorbarXPosLesion(temp)-1)),colorbarYPosLesion,:) = 0;

            else
    
                % map color bar to [0:1] ([0%:100%])
                colorbarXPercentLesion = -100 + 200*(colorbarXPosLesion-min(colorbarXPosLesion))/(max(colorbarXPosLesion)-min(colorbarXPosLesion));
                colorbarXRangeLesion = 128.5 + 127.5 * colorbarXPercentLesion/100;
                colorbarXColorsLesion = round(256*(interp1(1:256,colorMap,colorbarXRangeLesion)));
                % paint the colorbar
                for iYPos = 1: numel(colorbarYPosLesion)
                    canvas(size(canvas,1)-colorbarXPosLesion,colorbarYPosLesion(iYPos),:) = colorbarXColorsLesion;
                end
                % percent scale
                colorbarXPercentLabelsLesion = -100:10:100;
                colorbarXPercentLabelsPosLesion = colorbarXPosLesion(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), colorbarXPercentLabelsLesion));
                % add strokes to scale
                canvas( colorbarXPercentLabelsPosLesion, (colorbarYPosLesion(end)+1):(colorbarYPosLesion(end)+5), :) = 255;
        
                % remove color from scale of less than abs(params.cutoffDifferenceImage)
                lowerCutoff = -params.cutoffDifferenceImage * 100;
                upperCutoff = params.cutoffDifferenceImage * 100;
                mask = [lowerCutoff, upperCutoff];
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), mask, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(temp(1))):(colorbarXPosLesion(temp(2))-1)),colorbarYPosLesion,:) = 0;
                % remove color from scale above max percentage
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), max(uniquegbcminusgac)*100, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(temp):colorbarXPosLesion(end)-1)),colorbarYPosLesion,:) = 0;
                % remove color from scale below min percentage
                temp = cell2mat(arrayfun(@(x) find(min(abs(colorbarXPercentLesion-x)) == abs(colorbarXPercentLesion-x)), min(uniquegbcminusgac)*100, 'UniformOutput', false));
                canvas(size(canvas,1)-((colorbarXPosLesion(1):colorbarXPosLesion(temp)-1)),colorbarYPosLesion,:) = 0;
            end



            fig = figure;
            imshow(canvas)%, 'Border', 'tight')
            % add orientation and slice labels
            for r = 1:size(layout.data,1)
                for c = 1: size(layout.data,2)
                    text(layout.position{r,c}(1) + layout.position{r,c}(3) -20, layout.position{r,c}(2) + 20, layout.sliceLabel{r,c}, 'FontSize', 16, 'Color', [1 1 1], 'Units', 'data', 'HorizontalAlignment','right', 'VerticalAlignment', 'top')%, 'BackgroundColor', [0, 0, 0])
                    text(layout.position{r,c}(1) + 20, layout.position{r,c}(2) + 20, layout.orientationLabel{r,c}, 'Color', [1 1 1], 'FontSize', 16, 'Units', 'data', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')%, 'BackgroundColor', [0, 0, 0])
                    % add region labels on enlarged plots
                    if layout.doEnlarge(r,c) && ~params.RegionDetailsMonochrome
                        lineText = 1;
                        maxTextLength = max(cellfun(@numel,textstr));
                        % which labels are visible on this slice?
                        RoiDetailUsed = layout.RoiDetailUsed{r,c};
                        % display if either the left or the right region is visible
                        RoiDetailUsed = RoiDetailUsed(1:numel(RoiDetailUsed)/2) | RoiDetailUsed((numel(RoiDetailUsed)/2)+1:end);
                        for iROI = find(RoiDetailUsed)
                            while numel(textstr{iROI}) < maxTextLength
                                textstr{iROI} = [textstr{iROI}, ' '];
                            end
                            % background color is color of the ROI outline
                            text(layout.position{r,c}(1) + 20, layout.position{r,c}(2) + 60 + lineText * 24, textstr{iROI}, 'Color', [1,1,1], 'FontName', 'FixedWidth', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none', 'BackgroundColor', textcol{iROI})
                            lineText = lineText + 1;
                        end
                    end
                end
            end

            % add percent labels to colorbar
            for iLabel = 1:2:numel(colorbarXPercentLabelsLesion)
                text(colorbarYPosLesion(end)+10, size(canvas,1) - colorbarXPercentLabelsPosLesion(iLabel),  sprintf('%i %%', colorbarXPercentLabelsLesion(iLabel)), 'FontSize', 12, 'Color', [1,1,1])
            end
            if iDeficit == 5  
            % add n labels to colorbar
                tempMax = max(unique(gac));
                for iLabel = 1: 2: numel(colorbarXSubjectLabelsLesion)
                    if iLabel <= round(max(uniquegbcminusgac)*numel(lesionFilenames))+1 % +1 because starting at label "0"
                        text(colorbarYPosLesion(1)-10, size(canvas,1) - colorbarXSubjectPercentLabelsPosLesion(iLabel),  sprintf('%i', colorbarXSubjectLabelsLesion(iLabel)), 'FontSize', 12, 'Color', [1,1,1], 'HorizontalAlignment', 'right')
                    else
                        text(colorbarYPosLesion(1)-10, size(canvas,1) - colorbarXSubjectPercentLabelsPosLesion(iLabel),  sprintf('%i', colorbarXSubjectLabelsLesion(iLabel)), 'FontSize', 12, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment', 'right')
                    end
                
                end
            end
            if iDeficit == 5
                saveName = [figFolder, filesep, sprintf('LesionOverlap_combined_visThresh_%0.2f', params.cutoffDifferenceImage)];
            else
                saveName = [figFolder, filesep, sprintf('%s_threshold_%0.1f_combined_visThresh_%0.2f',params.criterionVariable, params.criterionThreshold, params.cutoffDifferenceImage)];
            end
            if params.doThalamusMasking
                saveName = [saveName, '_masked'];
            end
            if params.doFlipping
                saveName = [saveName, '_flipped'];
            end
            if params.replaceColorMap
                saveName = [saveName, sprintf('_ColorSaturationAt%0.2f', params.colorSaturationAt)];
            end
            if params.RegionDetailsMonochrome
                saveName = [saveName, '_ROIDetailsMonoChrome'];
            end
            if params.showOnlyRegionsDetails
                saveName = [saveName, '_ROIDetailsOnly'];
            end
            %saveName = [saveName, '_', sprintf('_%s', string(datetime('now','format','yyyy_MM_dd_HH_mm')))];
            axis equal
            axis tight
            axis off
            truesize

            savefig(fig, [saveName, '.fig']);
  
        end
    end
end
