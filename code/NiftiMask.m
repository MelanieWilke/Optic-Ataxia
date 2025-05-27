classdef NiftiMask < Nifti
    %NiftiMask.m class to represent a mask stored as nifti file
    %   Detailed explanation goes here

    properties
        Level
        Side
        ROIName
        Pattern = Pattern.Border;
        ColorFront = [255, 255, 255];
        ColorBack = []; 
    end

    methods
        function obj = NiftiMask(filename, name, level, pattern, colorFront, colorBack)
            
            % pre initialization


           
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
            args{3} = true;
            % superclass constructor
            obj = obj@Nifti(args{:});
            
            % post initialization 
            if nargin >= 5
                obj.ColorBack = colorBack;
            else
                obj.ColorBack = [];
            end
            if nargin >= 4
                obj.ColorFront = colorFront; 
            else
                obj.ColorFront = [255 255 255];
            end
            if nargin >= 3
                obj.Pattern = pattern;
            else
                obj.Pattern = Pattern.Border;  
            end
            if nargin >= 3
                obj.Level = level;
            else
                obj.Level = nan;
            end
            temp = regexp(obj.Name,'_','split');
            if numel(temp)==2
                obj.Side = temp{1};
                obj.ROIName = temp{2};
            elseif numel(temp)==3
                obj.Side = temp{1};
                obj.ROIName = strcat(temp{2}, '_', temp{3});
            else
                obj.Side = '';
                obj.ROIName = '';
            end
            
        end

        function img = slice(obj, direction, tal)
            img = slice@Nifti(obj, direction, tal);
            img.data = img.data(:,:,1);
            img.pattern = obj.Pattern;
            img.colorFront = obj.ColorFront;
            img.colorBack = obj.ColorBack;
            img.name = obj.Name;
            img.roiName = obj.ROIName;
            img.side = obj.Side; 
        end
    end
end