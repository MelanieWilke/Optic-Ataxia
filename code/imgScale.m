function imgS = imgScale(img, factor)
%imgScale enlarge the image by an integer factor
%   Detailed explanation goes here
if nargin < 2
    factor = 2;
end
if ~(round(factor)==factor) || factor < 1
    error('imgScale: factor needs to be an integer larger or equal to 1.')
end
imgS = img;
imgS.data = cell2mat(arrayfun(@(x) repmat(x,factor,factor),imgS.data,'UniformOutput',false));
% update T
for d = int32(Direction.X):int32(Direction.Z)
    if any(img.rc == d)
        % adjust midpoint of first voxel
        imgS.T(4,d) = imgS.T(4,d) + (0.5/factor - 0.5) * imgS.T(d,d);
        % stepsize is reduced by factor
        imgS.T(d,d) = imgS.T(d,d)/factor;
    end
end
