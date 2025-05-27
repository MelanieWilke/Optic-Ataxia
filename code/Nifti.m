classdef Nifti
    %Nifti class to represent a Nifti dataset
    %   Detailed explanation goes here

    properties
        Filename
        Name 
        Header
        Data
        isMask
    end

    methods
        function obj = Nifti(filename, name, isMask)
            %Nifti Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 3
                obj.isMask = isMask;
            else
                obj.isMask = false;
            end
            if nargin >= 2
                obj.Name = name;
            else
                obj.Name = '';
            end
            if nargin >= 1 && isfile(filename)
                obj.Filename = filename; 
                obj.Header = niftiinfo(filename);
                if obj.isMask
                    obj.Data = logical(niftiread(filename));
                else
                    obj.Data = niftiread(filename);
                end
            else
                obj.Filename = '';
                obj.Header = struct();
                obj.Data = [];
            end
        end

        function newobj = indexSelect(obj,indexRange)
            X = 1;
            Y = 2;
            Z = 3;
            directions = Direction(X:Z);
            for i = X:Z
                indexRange(i,:) = sort(indexRange(i,:));
                indexRange(i,1) = floor(indexRange(i,1));
                indexRange(i,2) = ceil(indexRange(i,2));
                if any(indexRange(i,:) < 1) || any(indexRange(i,:) > size(obj.Data, directions(i)))
                    error('Nifti.select): direction %s (mni [%i, %i], index [%i, %i]) is outside of the valid range [1, %i]', string(directions(i)), mniRange(i,1), mniRange(i,2), indexRange(i,1), indexRange(i,2), size(obj.Data, directions(i)));
                end
            end
            newobj = obj;
            newobj.Data = newobj.Data(indexRange(X,1):indexRange(X,2),indexRange(Y,1):indexRange(Y,2),indexRange(Z,1):indexRange(Z,2));
            for i = X:Z
                newobj.Header.Transform.T(4,i) = obj.index2mni(directions(i), indexRange(i,1));
            end
            newobj.Header.ImageSize = size(newobj.Data);    
        end
        
        function newobj = mniSelect(obj,mniRange)
            X = 1;
            Y = 2;
            Z = 3;
            directions = Direction(X:Z);
            indexRange = nan(size(mniRange));
            for i = X:Z
                indexRange(i,:) = sort(obj.mni2index(directions(i),mniRange(i,:)));
                indexRange(i,1) = floor(indexRange(i,1));
                indexRange(i,2) = ceil(indexRange(i,2));
                if any(indexRange(i,:) < 1) || any(indexRange(i,:) > size(obj.Data, directions(i)))
                    error('Nifti.select): direction %s (mni [%i, %i], index [%i, %i]) is outside of the valid range [1, %i]', string(directions(i)), mniRange(i,1), mniRange(i,2), indexRange(i,1), indexRange(i,2), size(obj.Data, directions(i)));
                end
            end
            newobj = obj;
            newobj.Data = newobj.Data(indexRange(X,1):indexRange(X,2),indexRange(Y,1):indexRange(Y,2),indexRange(Z,1):indexRange(Z,2));
            for i = X:Z
                newobj.Header.Transform.T(4,i) = obj.index2mni(directions(i), indexRange(i,1));
            end
            newobj.Header.ImageSize = size(newobj.Data);    
        end

        function newobj = enlarge(obj)
            newobj = obj.indexSelect([size(obj.Data)'/4, size(obj.Data)'*3/4]);
        end

        function img = slice(obj, direction, mni)
            index = obj.mni2index(direction, mni);
            % check index is valid
            img.index = round(index);
            if img.index < 1 || img.index > size(obj.Data, direction)
                error('Nifti.slice(): %i is outside of the valid range [1, %i] (%s direction of dataMatrix)', index, size(obj.Data, direction), string(direction));
            end
            img.data = [];
            % the numbers in the transpose matrix give the origin of the
            % first voxel in x, y, z.
            temp = sign(obj.Header.Transform.T(4,1:3));
            temp(temp==0)=1; % force 0 to be +0
            origin = [0:2:4]+(temp+3)/2;
            originLabels = {'L', 'P', 'I';...
                'R', 'A', 'S'};
            thisOriginLabels = Origin(origin);

            switch direction
                case Direction.X
                    img.data = repmat(squeeze(obj.Data(img.index,:,:)),1,1,3);
                    img.rc = [Direction.Y, Direction.Z];
                case Direction.Y
                    img.data = repmat(squeeze(obj.Data(:,img.index,:)),1,1,3);
                    img.rc = [Direction.X, Direction.Z];
                case Direction.Z
                    img.data = repmat(squeeze(obj.Data(:,:,img.index)),1,1,3);
                    img.rc = [Direction.X, Direction.Y];
            end
            img.rcOriginLabel= thisOriginLabels(img.rc);
           
            img.xy = img.rc(2:-1:1);

            img.originLabel = sprintf('%s%s', string(img.rcOriginLabel(1)),string(img.rcOriginLabel(2)));
            img.mni = obj.index2mni(direction, img.index);
            img.label = sprintf('%s = %s', string(direction), string(img.mni));
            img.T = obj.Header.Transform.T; 
            
        end

        function m = indexRange(obj)
            m = ones(3,2);
            m(:,2) = size(obj.Data)';
        end
        
        function m = mniRange(obj)
            X = 1;
            Y = 2;
            Z = 3;
            directions = Direction(X:Z);
            m = obj.indexRange();
            for i = X:Z
                m(i,:) = obj.index2mni(directions(i),m(i,:));
            end
        end
    end


    methods 
        function index = mni2index(obj, direction, mni)
            index = ((mni - obj.Header.Transform.T(4,direction))/ obj.Header.Transform.T(direction,direction)) + 1;
        end

        function mni = index2mni(obj, direction, index)
            mni = ((index - 1) * obj.Header.Transform.T(direction,direction)) + obj.Header.Transform.T(4,direction);
        end
    end
end
