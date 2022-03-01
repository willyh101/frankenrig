function [score allDistances]= scoreDistanceRules(matrix,varsToPass)

sz = size(matrix);
Dist = varsToPass.Dist;

targetRange = varsToPass.targetRange;
targDist = mean(targetRange);
if targDist == inf
    targDist = 2000;
elseif targDist == -inf
    targDist = 0;
end
targetRangeVal = (targetRange(2)-targetRange(1))/2;

minViolations=0; tooCloseScore =0; distScoreAgg = 0; 
for i = 1:sz(2)
    thisEns = matrix(:,i);
    subDist = Dist(thisEns,thisEns);
    
    minViolations = minViolations + sum(subDist(:)<varsToPass.minSep);
    
    minDistance  = min(subDist(:));
    tooCloseScore = tooCloseScore + varsToPass.minSep / minDistance; %drops to low number with more spread apart
    
    maxDistance = max(subDist(:));
    
    if maxDistance > targetRange(2)
        distScore = 100+(maxDistance-targetRange(2))^2;
    elseif maxDistance < targetRange(1)
        distScore = 100+(maxDistance-targetRange(1))^2;
    else
        distScore = abs(targDist-maxDistance) / targetRangeVal *100;
    end
    
    distScoreAgg = distScoreAgg + distScore;
    
    allDistances(i) = maxDistance;
end


    
score = minViolations*1000 + tooCloseScore*10 + distScoreAgg/10;
    