function gform = FormFunc(thickness, sccalib, dugd, dvgd)
% FORMFUNC Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;

%% group number
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);

%% thickness boundaries
groupboundary = [];
for ii=1:ngroup
    % unit: mm
    tmp = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.Thickness.Text);
    groupboundary = [groupboundary, tmp];
end    

%% Kernel form function groups
gform = [];
for ii=1:ngroup
    tmp = sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.ObjectScatterFit;
    % unit: cm -> mm
    sigma1 = str2double(tmp.sigma1.Text) * 10;
    sigma2 = str2double(tmp.sigma2.Text) * 10;
    B = str2double(tmp.B.Text);
    grid2 = dugd.^2 +dvgd.^2;
    gform(:,:,ii) = exp( -0.5 * grid2 /(sigma1^2)  ) + B * exp( -0.5 * grid2 /(sigma2^2) );
end
