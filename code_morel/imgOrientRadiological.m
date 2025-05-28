function img = imgOrientRadiological(img)
%imgOrientRadiological: orient image to radiological convention
%   Detailed explanation goes here

orthogonalDirection = setdiff(Direction(1:3), img.rc);
switch orthogonalDirection
    case Direction.X
        % swap rows to Direction.Z
        if img.rc(1) == Direction.Y
            img = imgSwapAxis(img);
        end
        % flip top to top
        if img.rcOriginLabel(1) == Origin.I
            img = imgFlipVertical(img);
        end
        % flip ant to left
        if img.rcOriginLabel(2) == Origin.P
            img = imgFlipHorizontal(img);
        end
    case Direction.Y
        % swap rows to Direction.Z
        if img.rc(1) == Direction.X
            img = imgSwapAxis(img);
        end
        % flip top to top
        if img.rcOriginLabel(1) == Origin.I
            img = imgFlipVertical(img);
        end
        % flip right to left
        if img.rcOriginLabel(2) == Origin.L
            img = imgFlipHorizontal(img);
        end
    case Direction.Z
        % swap rows to Direction.Y
        if img.rc(1) == Direction.X
            img = imgSwapAxis(img);
        end
        % flip anterior to top
        if img.rcOriginLabel(1) == Origin.P
            img = imgFlipVertical(img);
        end
        % flip right to left
        if img.rcOriginLabel(2) == Origin.L
            img = imgFlipHorizontal(img);
        end
end


end