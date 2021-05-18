function [outputArg1,outputArg2] = AmplitudeFactor(I0, Ip, sccalib)
%AMPLITUDEFACTOR Summary of this function goes here
%   Detailed explanation goes here

%% group number
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);

% Amplitude factor groups
cfactor = [];

for ii=1:ngroup
    tmp = sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.ObjectScatterFit;
    % Amplitude Factor
    % unit: cm^2 - > mm^2
    A = str2double(tmp.A.Text) * 100;
    % unitless
    alpha = str2double(tmp.alpha.Text);
    beta = str2double(tmp.beta.Text);
    cfactor(:,:,ii) = A.* (Ip./(I0+eps)).^(alpha) .* ( log(I0./(Ip+eps)) ).^(beta);    
end

end

