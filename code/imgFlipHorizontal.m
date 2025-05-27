function imgnew = imgFlipHorizontal(img)
%imgFlipHorizontal Flip the columns of an image 
% which are displayed on the x axis
%   Detailed explanation goes here
    % some constants
    ROW = 1;
    COLUMN = 2;
    imgnew = img;
    % flip the columns of the data
    nd = ndims(img.data);
    switch nd
        case 2 % mask
            imgnew.data = img.data(:,end:-1:1);
        case 3 % brain with color channels
            imgnew.data = img.data(:,end:-1:1,:);
        otherwise
            error('imgFlipVartical: img.data must be 3D (brain) or 2D (mask)');
    end
    % exchange the originLabel of the rows, i.e.
    % either flip R <-> L, A <-> P, or S <-> I
    % origins are stored in a matrix (numbered enumeration)
    % {'L'}{'R'}{'P'}{'A'}{'I'}{'S'}
    %  odd  even odd  even odd  even
    if rem(img.rcOriginLabel(COLUMN),2) % odd
        imgnew.rcOriginLabel(COLUMN) = img.rcOriginLabel(COLUMN) + 1; % change to even
    else % even
        imgnew.rcOriginLabel(COLUMN) = img.rcOriginLabel(COLUMN) - 1; % change to odd
    end
    % update originLabel 
    imgnew.originLabel = sprintf('%s%s', string(imgnew.rcOriginLabel(ROW)),string(imgnew.rcOriginLabel(COLUMN)));
    % update the T matrix
    nrows = size(imgnew.data, COLUMN);
    imgnew.T(4,imgnew.rc(COLUMN)) = imgnew.T(4,imgnew.rc(COLUMN)) + imgnew.T(imgnew.rc(COLUMN),imgnew.rc(COLUMN)) * (nrows-1);
    imgnew.T(imgnew.rc(COLUMN),imgnew.rc(COLUMN)) = imgnew.T(imgnew.rc(COLUMN),imgnew.rc(COLUMN)) * (-1);

end