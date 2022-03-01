function score = prioritizeStimmable(matrix,varsToPass)

roiWeights = varsToPass.roiWeights;

cells = unique(matrix(:));
cellWeights = roiWeights(cells); 

score = nanmean(cellWeights)*numel(matrix)/10;

