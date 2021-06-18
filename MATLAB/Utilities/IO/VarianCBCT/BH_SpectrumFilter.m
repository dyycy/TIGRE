function [outputArg1,outputArg2] = BH_SpectrumFilter(inputArg1,inputArg2)
%BH_SPECTRUMFILTER Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;




% I/I_0 = exp(- length * miu): dependent on incident energy
for ii = 1:length(sl)
    atten_tab(ii,:) = exp(-sl(ii).* BHCalib.bowtie.ac');
end

%%%%%%%%%%%%%%%->>>>>>>>>>>>>>>>>>>>>>
% kVp
kVp = length(BHCalib.source.spec);
% miu
BHCalib.object.ac =BHCalib.object.ac(1:kVp);
%% bowtie-attenuated spectra look-up table
% specLUT(sampling length, energy bin)
BHCalib.bowtie.specLUT = atten_tab(:,1:kVp).*BHCalib.source.spec;
%%%%%%%%%%%%%%%%%%%->>>>>>>>>>>>>>>




end

