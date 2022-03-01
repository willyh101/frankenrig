function [newMatrix] = updateMatrixFromList(startMatrix,varsToPass)

sz = size(startMatrix);

E = randperm(sz(2),1);
C = randperm(sz(1),1);

theseC = startMatrix(:,E);
newC = theseC(C);
iterUpdate = 0;
while ismember(newC,theseC)
    newC = varsToPass.cellsToUse(randperm(numel(varsToPass.cellsToUse),1));
    iterUpdate=iterUpdate+1;
    if iterUpdate >1000
        disp('Error did not find cell')
        break
    end
end

newMatrix = startMatrix;
newMatrix(C,E)=newC;
