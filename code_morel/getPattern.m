function [imgPattern1, imgPattern2] = getPattern(img, pattern)
%getPattern: get patterned img
% imgPattern1: img where pattern == true
% imgPattern2: img where pattern == false

if ~ismatrix(img)
    error('imgOutline: img is not a matrix')
end

patch = getPatch(pattern, size(img));
imgPattern1 = zeros(size(img),'like', img);
imgPattern1(patch) = img(patch);
imgPattern2 = zeros(size(img),'like', img);
if pattern ~= Pattern.Border
    imgPattern2(~patch) = img(~patch);
end