function label = generateLabel(label, msk, doEdge)
%generateLabel generate a patterned colored patch to represent the ROI in the legend  
%   Detailed explanation goes here
%label = reshape(repmat([0 0 0],c*r,1),r,c,3);
if nargin < 3
    doEdge = true;
end
mask = true(size(label,1), size(label,2));

% Face
[roiPattern1, roiPattern2] = getPattern(mask, msk.pattern);
label(repmat(roiPattern1>0,[1,1,3])) = reshape(repmat(msk.colorFront,[numel(label(roiPattern1>0)),1]),[],1);
if ~isempty(msk.colorBack)
    label(repmat(roiPattern2>0,[1,1,3])) = reshape(repmat(msk.colorBack,[numel(label(roiPattern2>0)),1]),[],1);
end
% Edge
if doEdge || msk.pattern == Pattern.Border
    labelBorder = getOutline(mask);
    label(repmat(labelBorder>0,[1,1,3])) = reshape(repmat(msk.colorFront,[numel(label(labelBorder>0)),1]),[],1);
end
end