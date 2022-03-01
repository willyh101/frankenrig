function [score ens_osi_vals] = scoreEnsTuning(matrix,varsToPass)

sz = size(matrix);
meanOriVals = varsToPass.meanOriVals;
calc_osi = @(po, oo) ((po - oo)/(po + oo));

clear meanTC
for i =1:sz(2)
    thisEns = matrix(:,i);
    thisOri = meanOriVals(thisEns,2:end);
    
    meanTC(i,:) = mean(thisOri);
end
meanTC = bsxfun(@minus, meanTC,min(meanTC,[],2));

[~, prefs] = max(meanTC, [], 2);

sz2=size(meanOriVals);
if sz2(2)==13;
    orthos1 = prefs-3;
    orthos1(orthos1<1)=orthos1(orthos1<1)+12;
    
    orthos2 = prefs+3;
    orthos2(orthos2>12)=orthos2(orthos2>12)-12;
elseif sz2(2)==9
    orthos1 = prefs-2;
    orthos1(orthos1<1)=orthos1(orthos1<1)+8;
    
    orthos2 = prefs+2;
    orthos2(orthos2>8)=orthos2(orthos2>8)-8;
else
    disp('ERROR UNEXPECTED SIZE')
end

orthos = [orthos1 orthos2];
clear osi_vals
for c=1:numel(prefs)
    po = prefs(c);
    oo = orthos(c,:);
    po_vals = meanTC(c, po);
    oo_vals = mean(meanTC(c,oo));
    ens_osi_vals(c) = calc_osi(po_vals, oo_vals);
end



targetRange = varsToPass.ensOSIRange;
targEnsOSI = mean(targetRange);
if targEnsOSI == inf
    targEnsOSI = 1;
elseif targEnsOSI == -inf
    targEnsOSI = 0;
end
targetRangeVal = (targetRange(2)-targetRange(1));

ensOSIScoreAgg =0;
for i = 1:sz(2)
    ensOSIVal = ens_osi_vals(i);
    if ensOSIVal > targetRange(2)
        ensOSIScore = 100+(abs(ensOSIVal-targetRange(2))/targetRangeVal*5)^2;
    elseif ensOSIVal < targetRange(1)
        ensOSIScore = 100+(abs(ensOSIVal-targetRange(1))/targetRangeVal*5)^2;
    else
        ensOSIScore = abs(targEnsOSI-ensOSIVal) / targetRangeVal *100;
    end
    ensOSIScoreAgg =ensOSIScoreAgg + ensOSIScore;
end

score = ensOSIScoreAgg; 
