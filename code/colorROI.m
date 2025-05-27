function imgR = colorROI(img, roi, doEdge)
%colorRoi apply patterened color to ROI
%   Detailed explanation goes here
if nargin < 3
    doEdge = true;
end
%some checks
if ~(all(img.T == roi.T,"all") && all(img.rc == roi.rc) && img.mni == roi.mni && numel(img.data) == numel(roi.data)*3) 
    error('colorRoi: mismatch between img and roi')
end

% some more to do

imgR = img;
% Face
[roiPattern1, roiPattern2] = getPattern(roi.data, roi.pattern);
imgR.data(repmat(roiPattern1>0,[1,1,3])) = reshape(repmat(roi.colorFront,[numel(img.data(roiPattern1>0)),1]),[],1);
if ~isempty(roi.colorBack)
    imgR.data(repmat(roiPattern2>0,[1,1,3])) = reshape(repmat(roi.colorBack,[numel(img.data(roiPattern2>0)),1]),[],1);
end
% Edge
if doEdge || roi.pattern == Pattern.Border
    roiBorder = getOutline(roi.data);
    imgR.data(repmat(roiBorder>0,[1,1,3])) = reshape(repmat(roi.colorFront,[numel(img.data(roiBorder>0)),1]),[],1);
end
end