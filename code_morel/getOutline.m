function imgOutline = getOutline(img)
%imgOutline get the outline of an image
%   i.e. all voxels that are at the border
% including voxels at the border of the image
if ~ismatrix(img)
    error('imgOutline: img is not a matrix')
end
% pad the image with zeros for border detection
imgMask = false(size(img,1)+2, size(img,2)+2);
imgMask(2:end-1,2:end-1) = (img ~= 0);
% use 8 border
shift = [repmat((-1:1)',3,1), sort(repmat((-1:1)',3,1))];
shift = shift(~all(shift==0,2),:); % remove [0, 0]
% detect if voxel has any false neighbour
isBorder = false(size(img));
for i = 1:size(shift,1) 
    isBorder = isBorder | ~imgMask(2+shift(i,1):end-1+shift(i,1),2+shift(i,2):end-1+(shift(i,2)));
end
% return border voxel
imgOutline = zeros(size(img),"like", img);
imgOutline(isBorder) = img(isBorder);