% defined beforehand
image = imread('./Datasets/corridor/001.png');
imageSize = size(image);
imGrid = struct('numCells', [6;5], 'skipOuter', false);

% create in 'createSensors.m'
imGrid.imPatches        = cell(imGrid.numCells(1), imGrid.numCells(2));
imGrid.imPatcheIndices  = cell(imGrid.numCells(1), imGrid.numCells(2));

patchHeight  = ceil(imageSize(1) / imGrid.numCells(1));
patchWidth   = ceil(imageSize(2) / imGrid.numCells(2));

for row = 1:imGrid.numCells(1)
    for col = 1:imGrid.numCells(2)
        
        x1 = max((col-1) * patchWidth,  1);
        x2 = min(x1 + patchWidth, imageSize(2));
        y1 = max((row-1) * patchHeight, 1);
        y2 = min(y1 + patchHeight, imageSize(1));
        
        imGrid.imPatcheIndices{row, col} = struct('x1', x1, 'x2', x2, ...
                                                  'y1', y1, 'y2', y2);
                                              
        imGrid.imPatches{row, col} = image(y1:y2, x1:x2);
    end
end