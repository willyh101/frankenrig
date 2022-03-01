function [score meanOSIs] = scoreMeanOSI(matrix,varsToPass)

sz = size(matrix);
OSIs = varsToPass.osi_vals;

clear meanOSIs;
for i=1:sz(2)
    thisEns = matrix(:,i);
    meanOSIs(i) = mean(OSIs(thisEns));
end

targetRange = varsToPass.meanOSIrange;
targEnsOSI = mean(targetRange);
if targEnsOSI == inf
    targEnsOSI = 1;
elseif targEnsOSI == -inf
    targEnsOSI = 0;
end
targetRangeVal = (targetRange(2)-targetRange(1));

ensOSIScoreAgg =0;
for i = 1:sz(2)
    ensOSIVal = meanOSIs(i);
    if ensOSIVal > targetRange(2)
        ensOSIScore = (abs(ensOSIVal-targetRange(2))/targetRangeVal*5)^2;
    elseif ensOSIVal < targetRange(1)
        ensOSIScore = (abs(ensOSIVal-targetRange(1))/targetRangeVal*5)^2;
    else
        ensOSIScore = 0; %abs(targEnsOSI-ensOSIVal) / targetRangeVal *100;
    end
    ensOSIScoreAgg =ensOSIScoreAgg + ensOSIScore;
end

score = ensOSIScoreAgg; 
