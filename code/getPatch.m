function patch = getPatch(pattern, sz)
%getPatch generate a pixel pattern
%   patterns are enumerated in Pattern.m

switch pattern
    case Pattern.Filled 
        patch = [1];
    case Pattern.Border
        patch = [0];
    case Pattern.Horizontal
        patch = [1;0];
    case Pattern.Vertical
        patch = [1 0];
    case Pattern.ObliqueLeft
        patch = [1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1];
    case Pattern.ObliqueRight
        patch = [0 0 1 1; 0 1 1 0; 1 1 0 0; 1 0 0 1];
    case Pattern.Pepper
        patch = [1 0; 0 1];
    case Pattern.Cross
        patch = [0 1 0; 1 1 1; 0 1 0];
end
r = ceil(sz(1)/size(patch,1));
c = ceil(sz(2)/size(patch,2));
patch = repmat(patch,r,c);
patch = logical(patch(1:sz(1), 1:sz(2)));
end