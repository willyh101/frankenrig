function makeVisStim()

global ExpStruct Exp_Defaults

visStartTime = 1000;

ExpStruct.triggerSI5=zeros(size(ExpStruct.triggerSI5))...
    + makepulseoutputs(1,1,25,1,1,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);

updateAOaxes();

blankOutput=zeros(size(ExpStruct.StimLaserEOM));

ExpStruct.triggerPuffer=blankOutput;
ExpStruct.nextholoTrigger = blankOutput;
ExpStruct.StimLaserEOM=blankOutput;
ExpStruct.nextsequenceTrigger = blankOutput;
ExpStruct.nextsequenceTrigger=zeros(size(ExpStruct.nextsequenceTrigger))...
            + makepulseoutputs(visStartTime,1,10,1,50,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);

updateAOaxes();
saveoutput('Out 0. 0mW. No Pulses');
disp('Saved Control Output');

