
function [finalMatrix score] = genericGDparfor(startMatrix,iterLim,updateFun,scoreFun,varsToPass,verbose)
%% Generic Gradient Descent

% verbose =1;
parCount = 5;

matrix = startMatrix;
score = scoreFun(startMatrix,varsToPass);

iterCount = 0;
while iterCount<iterLim
    iterCount = iterCount+ceil(parCount/2);
    [oldscore] = scoreFun(matrix,varsToPass);
    
    testMatrices = zeros([size(matrix) 5]);
    testScores = nan([1 5]);
    parfor i= 1:parCount
    newMatrix = updateFun(matrix,varsToPass);
    newScore = scoreFun(newMatrix,varsToPass); 
    
    testMatrices(:,:,i) = newMatrix;
    testScores(i) = newScore;
  
    end
    
    if any(testScores < oldscore)
        [newScore idx] = min(testScores);
        
        matrix = testMatrices(:,:,idx);
        score = newScore;
%         iterCount = 0;
        if verbose
            disp(num2str(newScore))
        end
        
    end
    
    if mod(iterCount,1001)==1000
        fprintf('.')
    end
end
finalMatrix = matrix; 
disp('done')
