% PlotMorelAtlas.m
% matlab script to generate a figure with axial and coronal slices of the
% Morel atlas.
% 
% run ProcessMorelAtlas.m first

% possible Layouts
asInFig5 = 1;
transposed = 2;

% choose layout
Layout = transposed;


% setting up more restricted fields of view
% magnification at the level of the thalamus
direction = [repmat(Direction.Y,1,5), repmat(Direction.Z,1,5)]; % Direction.Y, Direction.Z
pos = [-13, -16, -19, -22, -25, 2, 5, 8, 11, 14]; % MNI slices 
fovMNI = [-35, 35; -40, 10; -20, 30]; % field of view in MNI coordinates [xmin, xmax; ymin, yMax; zmin, zmax]

sz = size(mni.Data);
% setup for chosen layout
if Layout == asInFig5
    r = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
    c = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5];
    lblSize = [20, 40]; % size of the label patch
    lblFontSize = 12;
    % use this for larger fov as in figure  
    fovIndex = [round(sz(1)/4), round(sz(1)*3/4); round(sz(2)/4), round(sz(2)*3/4); round(sz(3)/4), round(sz(3)*3/4)];
    for d = 1:3
        for i = 1:2
            fovMNI(d,i) = mni.index2mni(Direction(d),fovIndex(d,i));
        end
    end
elseif Layout == transposed 
    c = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
    r = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5];
    lblSize = [24, 40]; % size of the label patch
    lblFontSize = 14;
else
    errore("unknown Layout")
end

% further settings to control the appearance of the figure
doRightOnly = true; % only right hemisphere of thalamus for sagittal slices
doExclusive = true; % exclusion of overlap between regions
doEdge = true; % plot a line around the ROI
doLeftEdgeOnly = true; %on the left thalamus show only outline of ROIs
Magnification = 2; % 1, 2, 3, 4 times magnification of images
doLevel = false; % plot ROIs of a level vs. plot Kumar 2022 ROIs
level = 3;


%% % deprecated code:
% % old layouts 
% switch direction
%     case Direction.Z
%         % display these axial slices
%         pos = [-12, -9, -6, -3, 0, 3, 6, 9, 12, 15, 18, 21];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%          % display these axial slices
%         pos = [-4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%         pos = [-13, -10, -7, -4, -1, 2, 5, 8, 11, 14, 17, 20];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%     case Direction.Y
%         pos = [-3, -6, -9, -12, -15, -18, -21, -24, -27, -30, -33, -36];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%         pos = [- 7, -9, -11, -13, -15, -17, -19, -21, -23, -25, -27, -29];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%          pos = [-1, -4, -7, -10, -13, -16, -19, -22, -25, -28, -31, -34];
%         % plot orientation
%         r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%         c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%     case Direction.X
%         if doRightOnly
%             pos = [ 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24];
%             % plot orientation
%             r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%             c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%         else
%             pos = [-20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 24];
%             % plot orientation
%             r = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
%             c = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
%         end
% end

%% choose regions
[~,index]=sort(order); % go through the ROIs in this order, choose ROI if in use

%% prepare slices
clear im
%im(max(r), max(c)) = struct();
clear msk;
rsize = zeros(max(r), max(c));
csize = zeros(max(r), max(c));
for i = 1:numel(pos)
    % mni
    % start with field of view
    thisfovMNI = fovMNI;
    % but use only a single slice in the direction you want to plot
    thisfovMNI(direction(i),:) = pos(i);
    % chose slab of data from the mni brain using mni coordinates from fov
    thisMNI = mni.mniSelect(thisfovMNI);
    % slice, orient to radiological convention and magnify
    im(r(i),c(i)) = imgScale(imgOrientRadiological(thisMNI.slice(direction(i), pos(i))), Magnification);
    tempSize = size(im(r(i),c(i)).data);
    rsize(r(i),c(i)) = tempSize(1);
    csize(r(i),c(i)) = tempSize(2);

    % ROIs
    if doLevel
        % ROIs by level
        for m = 1:numel(roiExclusiveByLevel{level})
            if doExclusive % non-overlapping ROIs
                % treatment as above for mni
                tmp = roiExclusiveByLevel{level}(m).mniSelect(thisfovMNI);
                msk(i,m) = imgScale(imgOrientRadiological(tmp.slice(direction(i), pos(i))), Magnification);
            else
                tmp = roiByLevel{level}(m).mniSelect(thisfovMNI);
                msk(i,m) = imgScale(imgOrientRadiological(tmp.slice(direction(i), pos(i))), Magnification);
            end
            if strcmp(msk(i,m).side, 'L') && doLeftEdgeOnly
                msk(i,m).pattern = Pattern.Border;
            end
            % do coloring of the masks only after all images have been reoriented,
            % because reorienting can interfere with the orientation of the
            % patterns which will then not match the pattern labels
            im(r(i),c(i)) = colorROI(im(r(i),c(i)),msk(i,m), doEdge);
        end
    else
        % ROIs by order and use (as in Kumar 2022)
        iAcc = 0;
        for m = 1:numel(roi)
            if use(index(m))== 1
                iAcc = iAcc +1;
                if doExclusive
                    tmp = roiExclusive(index(m)).mniSelect(thisfovMNI);
                    msk(i,iAcc) = imgScale(imgOrientRadiological(tmp.slice(direction(i), pos(i))), Magnification);
                else
                    tmp = roi(index(m)).mniSelect(thisfovMNI);
                    msk(i,iAcc) = imgScale(imgOrientRadiological(tmp.slice(direction(i), pos(i))), Magnification);
                end
                if strcmp(msk(i,iAcc).side, 'L') && doLeftEdgeOnly
                    msk(i,iAcc).pattern = Pattern.Border;
                end
                % do coloring of the masks only after all images have been reoriented,
                % because reorienting can interfere with the orientation of the
                % patterns which will then not match the pattern labels
                 im(r(i),c(i)) = colorROI(im(r(i),c(i)),msk(i,iAcc), doEdge);
            end
        end
    end
end

%% canvas
% cave! will not work for empty cells

% calculate positions for images (slices) and labels
% images are oriented in a grid
% for every row and column get the max image size
max_rsize = max(rsize,[],2);
max_csize = max(csize);

patch_r_end = cumsum(max_rsize);
patch_r_start = [1; patch_r_end(1:end-1)+1];
patch_c_end = cumsum(max_csize);
patch_c_start = [1, patch_c_end(1:end-1)+1];

% if some images are smaller: center them by calculating an offset 
r_offset = table2array(varfun(@(x) floor((max_rsize-x)/2),table(rsize)));
c_offset = table2array(rowfun(@(x) floor((max_csize-x)/2),table(csize)));

% repmat starts and ends to make them r
patch_r_start = repmat(patch_r_start, 1, max(c));
patch_r_end = repmat(patch_r_end, 1, max(c));
patch_c_start = repmat(patch_c_start, max(r), 1);
patch_c_end = repmat(patch_c_end, max(r), 1);
% calculate image positions
image_r_start = patch_r_start + r_offset;
image_r_end = image_r_start + rsize - 1;
image_c_start = patch_c_start + c_offset;
image_c_end = image_c_start + csize - 1;


% startIndices = [1,1; round(imsize(1:2)/4); round(imsize(1:2)/3)+[10,0]];
lblRStart = lblSize(1);
lblRStep = lblSize(1) + 7;
%lblCStart = cumsum_max_csize(end) + 20;
lblCStart = patch_c_end(end) + 20;
txtRStart = lblRStart + lblSize(1);
txtCStart = lblCStart + lblSize(2) + 10;
lblSpace = lblRStart+((size(msk,2)/2)+1)*lblRStep;
% prepare canvas
canvas = 127*ones(max(max(patch_r_end,[],'all'), lblSpace), max(patch_c_end,[],'all') + 150, size(im(1).data, 3), 'like', im(1).data);
% insert images
imgpos = [];
for i = 1:numel(pos)
    canvas(image_r_start(i): image_r_end(i), image_c_start(i): image_c_end(i),:) = im(i).data;
end

% labels and text (use only right ROIs)
m = 0;
for n = 1:size(msk,2)
    if msk(1,n).side == 'R'
        m = m+1;
        canvas(lblRStart+1+(m-1)*lblRStep: lblRStart+lblSize(1)+(m-1)*lblRStep, lblCStart: lblCStart + lblSize(2) - 1, :) = generateLabel(canvas(lblRStart+1+(m-1)*lblRStep: lblRStart+lblSize(1)+(m-1)*lblRStep, lblCStart: lblCStart + lblSize(2) - 1, :), msk(1,n), doEdge);
    end
end
% display
fig = figure;
image(canvas)
% add image orientation label and slice number
for i = 1:numel(pos)
    t1(i) = text(image_c_start(i) + 20, image_r_start(i) + 20, im(i).originLabel(end:-1:1), Color=[1,1,1], FontSize=16, VerticalAlignment="top");
    t2(i) = text(image_c_end(i) - 20, image_r_start(i) + 20, lower(im(i).label), Color=[1,1,1], FontSize=16, HorizontalAlignment="right", VerticalAlignment="top");
end
% add ROI labels
m= 0;
for n = 1:size(msk,2)
    if msk(1,n).side == 'R'
        m = m+1;
        t(m) = text(txtCStart,txtRStart + (m-1)*lblRStep, regexprep(msk(1,n).roiName,'_',' '), FontSize = lblFontSize, Color = [0,0,0], VerticalAlignment="bottom");
    end
end

ax = gca;
axis equal
axis tight
axis off
truesize

