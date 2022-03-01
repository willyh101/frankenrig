function [score, percentCoTuned] = scoreCoTuned(matrix,varsToPass)

prefDirs = varsToPass.POs;
oris = varsToPass.oris;

if max(oris)==330
    prefOris = prefDirs-6;
    prefOris(prefOris<2) = prefOris(prefOris<2)+6;
elseif max(oris)==315
    prefOris = prefDirs-4;
    prefOris(prefOris<2) = prefOris(prefOris<2)+4;
end

sz = size(matrix);

violations = 0;
for i=1:sz(2)
    thisEns = matrix(:,i);
    ensTunings = prefOris(thisEns);
    uTunings = unique(ensTunings);
    [~,~,ic] = unique(ensTunings);
    a_counts = accumarray(ic,1);
%     value_counts = [C', a_counts];
    [count,~] = max(a_counts);
%     tune = C(idx);
    percentUntuned = 1-count/sz(1);
    violations = violations + percentUntuned;
    percentCoTuned(i) = count/sz(1);
end

score = violations*1000;