
function out = outputPatternTranslator(ExpStruct,uniqueStims)

OutputStims=[];OutputNames = [];OutputOrder=[];OutputPatterns=[];OutputStimID=[];
for i =1:numel(uniqueStims)
    S=uniqueStims(i);
    
%     t = find(stimID==S,1);
    
    thisOutput = ExpStruct.stimlog{S}{1};
    
    isFound=0;
    c=0;
    while ~isFound
        c=c+1;
        
        if c>numel(ExpStruct.output_patterns)
            isFound = 1;
            disp('Output Not Found')
            c=1;
        else
            testOutput = ExpStruct.output_patterns{c};
            
            if numel(testOutput) == numel(thisOutput)
                isFound=  all(thisOutput(:)==testOutput(:));
                
            end
        end
        
    end
    
    str =strsplit(ExpStruct.output_names{c}(4:end),'.');
    
    
    OutputPatterns{i}= ExpStruct.output_patterns{c};
    OutputStims(i)=c;
    OutputNames{i}=ExpStruct.output_names{c};
    OutputStimID(i)=S;
    try
    OutputOrder(i) = str2num(str{1});
    catch
        OutputOrder(i) = 0;
        disp('Output Order Error, Likely because of non conventional Name')
    end

end

out.OutputStims     = OutputStims;      %number in stimLog
out.OutputNames     = OutputNames;      %name as written in output selector
out.OutputOrder     = OutputOrder;      %order to sort by to get back into stimparams order
out.OutputPatterns  = OutputPatterns;   %full pattern of outputs
out.OutputStimID    = OutputStimID;     %number in stimID