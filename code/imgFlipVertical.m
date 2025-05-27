function imgnew = imgFlipVertical(img)
%imgFlipVertical Flip the rows of an image 
% which are displayed on the y axis
%   Detailed explanation goes here
    % some constants
    ROW = 1;
    COLUMN = 2;
    imgnew = img;
    % flip the rows of the data
    nd = ndims(img.data);
    switch nd
        case 2 % mask
            imgnew.data = img.data(end:-1:1,:);
        case 3 % brain with color channels
            imgnew.data = img.data(end:-1:1,:,:);
        otherwise
            error('imgFlipVartical: img.data must be 3D (brain) or 2D (mask)');
    end
    % exchange the originLabel of the rows, i.e.
    % either flip R <-> L, A <-> P, or S <-> I
    % origins are stored in a matrix (numbered enumeration)
    % {'L'}{'R'}{'P'}{'A'}{'I'}{'S'}
    %  odd  even odd  even odd  even
    if rem(img.rcOriginLabel(ROW),2) % odd
        imgnew.rcOriginLabel(ROW) = img.rcOriginLabel(ROW) + 1; % change to even
    else % even
        imgnew.rcOriginLabel(ROW) = img.rcOriginLabel(ROW) - 1; % change to odd
    end
    % update originLabel 
    imgnew.originLabel = sprintf('%s%s', string(imgnew.rcOriginLabel(ROW)),string(imgnew.rcOriginLabel(COLUMN)));
    % update the T matrix
    nrows = size(imgnew.data, ROW);
    imgnew.T(4,imgnew.rc(ROW)) = imgnew.T(4,imgnew.rc(ROW)) + imgnew.T(imgnew.rc(ROW),imgnew.rc(ROW)) * (nrows-1);
    imgnew.T(imgnew.rc(ROW),imgnew.rc(ROW)) = imgnew.T(imgnew.rc(ROW),imgnew.rc(ROW)) * (-1);

end