function [score percentAboveThreshold]  = prioritizeVispResp(matrix,varsToPass)

p_vis_resp = varsToPass.p_vis_resp;

cells = unique(matrix(:));
cellpVis = max(log10(p_vis_resp(cells))-log10(varsToPass.pVisLim),0); 

score = nanmean(cellpVis)*numel(matrix)/5;

percentAboveThreshold = mean(cellpVis~=0);

