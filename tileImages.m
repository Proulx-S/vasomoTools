function tiledMatrix = tileImages(imToTile)
% TILEIMAGES - Tile multiple images into a single large matrix
%
% Input:
%   imToTile - Cell array of images to tile (all images must have same dimensions)
%
% Output:
%   tiledMatrix - Large matrix containing all images tiled in a grid
%                 Preserves all dimensions beyond the first two
%
% Example:
%   tiledMatrix = tileImages({img1, img2, img3, img4});
%   imagesc(tiledMatrix); axis image; colormap gray;

% Calculate optimal grid dimensions to be close to square
nImages = length(imToTile);
nCols = ceil(sqrt(nImages));
nRows = ceil(nImages / nCols);
sz = size(imToTile{1});

% Handle multi-dimensional images
if length(sz) == 2
    % 2D case - original behavior
    tiledMatrix = zeros(nRows * sz(1), nCols * sz(2));
    
    for i = 1:nImages
        % Calculate position in the large matrix
        rowStart = floor((i-1) / nCols) * sz(1) + 1;
        colStart = mod(i-1, nCols) * sz(2) + 1;
        
        % Place the image in the large matrix
        tiledMatrix(rowStart:rowStart+sz(1)-1, colStart:colStart+sz(2)-1) = imToTile{i};
    end
    
else
    % Multi-dimensional case (3D, 4D, etc.)
    % Create output with tiled first two dimensions, preserving others
    outputSize = [nRows * sz(1), nCols * sz(2), sz(3:end)];
    tiledMatrix = zeros(outputSize);
    
    for i = 1:nImages
        % Calculate position in the large matrix (first two dimensions only)
        rowStart = floor((i-1) / nCols) * sz(1) + 1;
        colStart = mod(i-1, nCols) * sz(2) + 1;
        
        % Create index arrays for all dimensions
        rowIdx = rowStart:rowStart+sz(1)-1;
        colIdx = colStart:colStart+sz(2)-1;
        
        % Handle different dimensionalities
        if length(sz) == 3
            % 3D case
            tiledMatrix(rowIdx, colIdx, :) = imToTile{i};
        elseif length(sz) == 4
            % 4D case
            tiledMatrix(rowIdx, colIdx, :, :) = imToTile{i};
        elseif length(sz) == 5
            % 5D case
            tiledMatrix(rowIdx, colIdx, :, :, :) = imToTile{i};
        else
            % General case for higher dimensions
            % Create a cell array of indices
            idx = cell(1, length(sz));
            idx{1} = rowIdx;
            idx{2} = colIdx;
            for d = 3:length(sz)
                idx{d} = ':';
            end
            tiledMatrix(idx{:}) = imToTile{i};
        end
    end
end

end 