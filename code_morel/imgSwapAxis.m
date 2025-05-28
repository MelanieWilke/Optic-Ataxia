function imgnew = imgSwapAxis(img)
%imgSwapAxis transpose the rows and columns of a matrix
%   Detailed explanation goes here
    ROW = 1;
    COLUMN = 2;
    COLOR = 3;
    imgnew = img;
    nd = ndims(img.data);
    switch nd
        case 2
            imgnew.data = permute(img.data,[COLUMN, ROW]);
        case 3
            imgnew.data = permute(img.data,[COLUMN, ROW, COLOR]);
        otherwise
            error('imgSwapAxis: img.data must be 3D (color) or 2D (mask)');
    end
    imgnew.rc = img.rc([COLUMN, ROW]);
    imgnew.rcOriginLabel = img.rcOriginLabel([COLUMN, ROW]);
    imgnew.xy = img.xy([COLUMN, ROW]);
    imgnew.originLabel = sprintf('%s%s', string(imgnew.rcOriginLabel(ROW)),string(imgnew.rcOriginLabel(COLUMN)));
end