function [score, minDistances] = scoreDistanceFromCenter(matrix,varsToPass)

sz = size(matrix);
DistFromCenter = varsToPass.DistFromCenter;
targetMinDist = varsToPass.minDistFromCenter;

minViolations = 0;
for i = 1:sz(2)
    thisEns = matrix(:,i);
    ensDist = DistFromCenter(thisEns');
    
    minViolations = minViolations + sum(ensDist<targetMinDist);
    minDistances(i) = min(ensDist);
end

score = minViolations*1000;