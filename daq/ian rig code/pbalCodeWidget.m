firstTimes = ExpStruct.bigListOfFirstStimTimes(:,1);

DAQEpoch = 4;
sw1 =ExpStruct.EpochEnterSweep{DAQEpoch};
try
    sw2 =ExpStruct.EpochEnterSweep{DAQEpoch+1};
catch
    sw2 = ExpStruct.sweep_counter-1;
end
swp = sw1:sw2;
powers = ExpStruct.outID(swp);

powerList = holoRequest.holoStimParams.powerList{1};
for i =1:numel(powers)
    if powers(i)>0
        powers(i)=powerList(powers(i));
    end
end


save('T:\holography\FrankenRig\OnlinePCurve\data.mat','powers','firstTimes', '-v7.3')