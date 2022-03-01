function [score] = patternScorer(pattern,target,maxCellPerBin)

score =0;

cellsPerBin = sum(pattern);

singles = sum(cellsPerBin==1);
doubles = sum(cellsPerBin==2);
overLim = sum(cellsPerBin>maxCellPerBin);

score = score + singles*25 + doubles*1000 + overLim*1000;

%% Try to optimize freq too

nCells = size(pattern,1);
minInter=zeros([1 nCells]);

for k=1:nCells
    mn = min(diff(find(pattern(k,:))));
    if isempty(mn)
        mn=0;
    end
    minInter(k) = mn;
end

temp = target;
temp(temp<=1)=NaN;

desiredInter = round(size(pattern,2)./(temp-1)-1);

interScore = nansum(abs(desiredInter./minInter-1));

closeInter = sum(minInter<3 & minInter~=0);
wayToCloseInter = sum(minInter<2 & minInter~=0);


score = score + interScore + closeInter*50 + wayToCloseInter*1000;