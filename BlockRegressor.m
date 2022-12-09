function [ActivationRegSelected,CapniaTrigSeries_delay] = BlockRegressor(Time, TriggerBlocks, Tau, Tau2, DelayTime, PlotOn, flag_no_SPION_decay, flag_no_hemodynamic_tau)
%This model assumes the HRF is two real poles with time constant Tau
%TriggerBlocks is 0 to 1 activation

F = @(T) T.*exp(-(T)/Tau); %making two functions that are the impulse response which describe the dynamics of the system
F2 = @(T) exp(-(T)/Tau2); %particle decay over time

RawBlocks = TriggerBlocks*2-1; %Makes it so the trigger is +1/-1 not 1/0

fs = 100; 
dt = 1/fs; 

TimeUpsample = 0;
BlocksUpsample = RawBlocks(1);
for i = 1:length(Time)-1
   Inds(i) = length(TimeUpsample); 
   TimeVec_i = (Time(i)+dt):dt:Time(i+1);
   if i == length(Time)-1 && Time(i+1) > TimeVec_i(end)
       TimeVec_i = (Time(i)+dt):dt:(Time(i+1)+dt);
   end  
   TimeUpsample = [TimeUpsample,TimeVec_i];
   LengthStep = length(TimeVec_i);
   BlocksUpsample = [BlocksUpsample,ones(1,LengthStep)*RawBlocks(i)];    
end
Inds = [Inds,length(BlocksUpsample)];

% time delay Blocks:
if DelayTime >= 0
    samples_to_delay = round(DelayTime/dt);
    BlocksUpsample = [BlocksUpsample(1)*ones(1,samples_to_delay),BlocksUpsample(1:end-samples_to_delay)];
else
    samples_to_delay = round(abs(DelayTime)/dt);
    BlocksUpsample = [BlocksUpsample(samples_to_delay+1:end),BlocksUpsample(end)*ones(1,samples_to_delay)];
end

% downsample delayed version for plotting: 
temp = interp1(TimeUpsample,BlocksUpsample,Time);
CapniaTrigSeries_delay = round((temp + 1)/2);


Timeend = 10*Tau;

TimeHRF = 0:dt:Timeend;
TimeDecay = TimeUpsample; 

if flag_no_SPION_decay
    Decay = ones(size(TimeDecay));
else
    Decay = F2(TimeDecay);
end

if flag_no_hemodynamic_tau
    HRF_TS = 1; 
else
    HRF_TS = F(TimeHRF);
end


BlockRegTmp = conv(BlocksUpsample,HRF_TS);
BlockReg = BlockRegTmp(1:end-length(HRF_TS)+1);

BlockReg = BlockReg/max(BlockReg);
BlockRegwDecay = BlockReg.*Decay;

ActivationRegSelected = BlockRegwDecay(Inds);

end