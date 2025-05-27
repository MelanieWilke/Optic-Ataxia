classdef NiftiWeight < Nifti
    %NiftiWeight: class to represent partially overlapping ROIs
    %   The Morel Atlas regions were found to be partially overlapping.
    %   Weights were used to represent shared voxels in a way that would
    %   not count shared voxels twice and that would keep sizes the same
    %   across the levels of the thalamic hierarchy.
    %   Overlapping voxels were assigned to all sharing regions with 
    %   reduced weight.
    %   e.g.:
    %       voxel in 1 ROI only: weight 1
    %       voxel shared between 2 ROIS: weight 0.5 in each ROI
    %       etc.
    %   Reduced weigths were propagated to higher levels of the Morel
    %   hierarchy
    %


    properties
        Level
        Side
        ROIName
    end

    methods
        function obj = NiftiWeight(filename, name, level)
            % pre initialization
            args{3} = false;
            if nargin >= 2
                args{2} = name;
            else
                args{2} = '';
            end
            if nargin >= 1
                args{1} = filename;
            else
                args{1} = '';
            end
            
            % superclass constructor
            obj = obj@Nifti( args{:});
            
            % post initialization
            if nargin == 3
                obj.Level = level;
            else
                obj.Level = nan;
            end
            temp = regexp(obj.Name,'_','split');
            if numel(temp)==2
                obj.Side = temp{1};
                obj.ROIName = temp{2};
            else
                obj.Side = '';
                obj.ROIName = '';
            end        
        end
    end
end