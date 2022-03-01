
function [finalMatrix score] = genericGD(startMatrix,iterLim,updateFun,scoreFun,varsToPass,verbose)
%% Generic Gradient Descent

% verbose =1;

matrix = startMatrix;
score = scoreFun(startMatrix,varsToPass);

iterCount = 0;
while iterCount<iterLim
    iterCount = iterCount+1;
    [oldscore] = scoreFun(matrix,varsToPass);
    
    newMatrix = updateFun(matrix,varsToPass);
    newScore = scoreFun(newMatrix,varsToPass); 
    
    if newScore < oldscore
        matrix = newMatrix;
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
