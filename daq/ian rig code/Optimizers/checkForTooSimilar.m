function score = checkForTooSimilar(matrix,varsToPass)

sz = size(matrix);

tooHighThresh = floor(sz(1)*varsToPass.percentSimilar/100);

c = 0; %count of duplicated cells
cTooHigh = 0; %count of too High Cells
for i = 1:sz(2)
    for k=1:sz(2)
        if i~=k
           numDuplicates = sum(ismember(matrix(:,i),matrix(:,k)));
           c = c + numDuplicates;
            if numDuplicates > tooHighThresh
            cTooHigh = cTooHigh+ numDuplicates-tooHighThresh;
            end
        end
    end
end
score = c + cTooHigh*1000;