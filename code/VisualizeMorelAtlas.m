% VisualizeMorelAtlas.m
% 
% matlab script to prepare visualization of the Morel atlas regions 
% on the MNI 152 T1 0.5mm template
%
% start from code folder
% 
% the visualization itself is done by a second script:
% plotMorelAtlas.m
%
% cf. README.md

clear;
close all;

thisFolder = pwd;
% check file was started correctly
lastFilesepPosPlus1 = regexp(thisFolder,'code$');
if ~isempty(lastFilesepPosPlus1)
    codeFolder = thisFolder;
    lastFilesepPos = lastFilesepPosPlus1 - 1; % pointer was at 'c', now is at filesep
    % generate name of dataFolder
    dataFolder = [codeFolder(1:lastFilesepPos), 'data'];
else
    error('please start VisualizeMorelAtlas from ...\Optic-Ataxia\code')
end

%% MNI template 0.5mm

MNITemplateFileName = [dataFolder, filesep, 'MNI152_T1_0.5mm.nii.gz'];

mni = Nifti(MNITemplateFileName, 'MNI152_T1_0.5mm');

%% Header transform
mni.Header.Transform.T
% 
%    -0.5000         0         0         0
%          0    0.5000         0         0
%          0         0    0.5000         0
%    90.0000 -126.0000  -72.0000    1.0000

% https://neuroimage.usc.edu/brainstorm/CoordinateSystems
% (X, Y, Z) axes are oriented towards (right, anterior, dorsal)
% first voxel in the matrix is right (X = 90), posterior, ventral

%% Image matrix size
size(mni.Data)
%    364   436   364
%    X     Y     Z

%% Morel atlas hierarchy
roiTable = readtable([dataFolder, filesep, 'Morel_hierarchical_2025_05.xlsx']);

roiNames = roiTable.regionNames;
roiLevels =roiTable.level;
use = roiTable.use;
order = roiTable.order;  
patt = roiTable.pattern;
r1 = roiTable.r1;
g1 = roiTable.g1;
b1 = roiTable.b1;
r2 = roiTable.b2;
g2 = roiTable.g2;
b2 = roiTable.b2;
for iROI = 1:numel(r1)
    if ~isnan(r1(iROI))
        colorFront{iROI} = round(255*[r1(iROI),g1(iROI),b1(iROI)]);
    else
        colorFront{iROI} = [];        
    end
    if ~isnan(r2(iROI))
        colorBack{iROI} = round(255*[r2(iROI),g2(iROI),b2(iROI)]);
    else
        colorBack{iROI} = [];        
    end
end

%% rois from Morel atlas
roi(numel(roiNames),1)=NiftiMask();
roiFileNames = ExtractAndFilterFiles([dataFolder, filesep, '\morel-atlas_2024_09_16'], 'nii.gz');
for iROI = 1:numel(roiNames)
    thisFileName = roiFileNames{cellfun(@(x) contains(x,strcat(roiNames{iROI},'.')),roiFileNames)};
    roi(iROI) = NiftiMask(thisFileName, roiNames{iROI}, roiLevels(iROI), Pattern(patt(iROI)), colorFront{iROI}, colorBack{iROI});
end


%% roiWeights from Morel atlas
roiWeight(numel(roiNames),1)=NiftiWeight();
roiWeightFileNames = ExtractAndFilterFiles([dataFolder, filesep, '\morel-atlas_2024_09_16_weights'], 'nii.gz');
for iROI = 1:numel(roiNames)
    thisFileName = roiWeightFileNames{cellfun(@(x) contains(x,strcat(roiNames{iROI},'_')),roiWeightFileNames)};
    roiWeight(iROI) = NiftiWeight(thisFileName, roiNames{iROI}, roiLevels(iROI));
end

%% rois from weights == 1
roiExclusive = roi;
for iROI = 1:numel(roi)
    roiExclusive(iROI).Data = roiWeight(iROI).Data == 1;
end

%% roi vs roiExclusive
for iROI = 1:numel(roi)
    fprintf('%i, %s: size roi: %i, size roiExclusive: %i\n', iROI, roiExclusive(iROI).Name, sum(roi(iROI).Data,'all'), sum(roiExclusive(iROI).Data, 'all'));               
end
%% sort by level
maxLevel = max(roiLevels); 
for iLevel = maxLevel:-1:1
    roiByLevel{iLevel} = roi(roiLevels==iLevel);
    roiExclusiveByLevel{iLevel} = roiExclusive(roiLevels==iLevel);
end

for iLevel = 1:maxLevel-1
    fprintf('Level %i\n',iLevel)
    for iROI = 1:numel(roiExclusiveByLevel{iLevel})
        fprintf('%s: %i\n', roiExclusiveByLevel{iLevel}(iROI).Name, sum(roiExclusiveByLevel{iLevel}(iROI).Data,"all"))
    end
end

for iLevel = 1:maxLevel-1
    fprintf('Level %i\n',iLevel)
    for iROI = 1:numel(roiExclusiveByLevel{iLevel})
        for iROIsBelow = 1: numel(roiExclusiveByLevel{iLevel+1})
            if sum(roiExclusiveByLevel{iLevel}(iROI).Data & roiExclusiveByLevel{iLevel+1}(iROIsBelow).Data, "all")
                fprintf('%s contains %s\n', roiExclusiveByLevel{iLevel}(iROI).Name, roiExclusiveByLevel{iLevel+1}(iROIsBelow).Name)
            end
        end
    end
end

for iLevel = 1:maxLevel
    M{iLevel} = nan(numel(roiByLevel{iLevel}),numel(roiByLevel{iLevel}));
    for iROI = 1: numel(roiByLevel{iLevel})
        fprintf('ROI: %s: voxel: %i, exclusive: %i (%.2f percent).\n',roiByLevel{iLevel}(iROI).Name, sum(roiByLevel{iLevel}(iROI).Data>0,'all'), sum(roiExclusiveByLevel{iLevel}(iROI).Data>0,'all'), round(100* sum(roiExclusiveByLevel{iLevel}(iROI).Data>0,'all')/sum(roiByLevel{iLevel}(iROI).Data>0,'all'),2))
        for jROI = 1:numel(roiExclusiveByLevel{iLevel})
            M{iLevel}(iROI,jROI) = sum(roiByLevel{iLevel}(iROI).Data & roiByLevel{iLevel}(jROI).Data, 'all');
        end
    end
    M{iLevel}   
end

%% build tree


% %% set color and pattern
% % level 2
% 
% ROINames2 = unique(arrayfun(@(x) x.ROIName, roiExclusiveByLevel{2}, 'UniformOutput', false));
% FaceColors2 = 255*[[1, 1, 0]; [1, 0, 0];  [0.2, 0.2, 1]; [0, 1, 0];];
% for iROI = 1:numel(roiExclusiveByLevel{2})
%     colorindex = cellfun(@(x) strcmp(x, roiExclusiveByLevel{2}(iROI).ROIName), ROINames2);
%     roiExclusiveByLevel{2}(iROI).ColorFront = FaceColors2(colorindex,:);
%     if strcmp(roiExclusiveByLevel{2}(iROI).Side, 'R')
%         roiExclusiveByLevel{2}(iROI).Pattern = Pattern.Filled;
%     else
%         roiExclusiveByLevel{2}(iROI).Pattern = Pattern.Border;
%     end
%     % pass colors to next level
%     for iROIsBelow = 1: numel(roiExclusiveByLevel{3})
%         if sum(roiExclusiveByLevel{2}(iROI).Data & roiExclusiveByLevel{3}(iROIsBelow).Data, "all")
%             fprintf('%s contains %s\n', roiExclusiveByLevel{2}(iROI).Name, roiExclusiveByLevel{3}(iROIsBelow).Name)
%             roiExclusiveByLevel{3}(iROIsBelow).ColorFront = roiExclusiveByLevel{2}(iROI).ColorFront;
%         end
%     end
% end
% 
% ROINames3 = unique(arrayfun(@(x) x.ROIName, roiExclusiveByLevel{3}, 'UniformOutput', false));


