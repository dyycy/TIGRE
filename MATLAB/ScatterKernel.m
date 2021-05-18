function [outputArg1,outputArg2] = ScatterCorrection(datafolder, Blk, BlkAirNorm, geo, ScanXML)
%SCATTERKERNEL Summary of this function goes here
%   Detailed explanation goes here

% unit: mm
offset=geo.offDetector;
% Center Kernel
offset = [0, 0];

% grid unit: mm
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% Downsampling
% unit: mm
% Downsampled to about 10 mm in axial direction
dus = downsample(us, 26);
% Downsampled to about 4 mm in axial direction
dvs = downsample(vs, 10);

step_du = mean(diff(dus));
step_dv = mean(diff(dvs));

%% Grid
% unit: mm
[ugd,vgd] = meshgrid(us,vs); %detector
[dugd, dvgd] = meshgrid(dus,dvs); %detector

%% Load Scatter Calibration
sccalib = ScCalibFromXML(datafolder);

%% Load Blk (comment out later)
%{
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);
% Detector point scatter correction
Blk = DetectorPointScatterCorrection(Blk, geo);
%}
sBlk = sum(Blk, 3);
sAirNorm = sum(BlkAirNorm);

%% n-thickness group number and boundaries
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);
nbounds = [];

for ii=1:ngroup
    % unit: mm
    tmp = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.Thickness.Text);
    nbounds = [nbounds, tmp];
end

%% Scatter Correction
for ii = 1:ScanXML






end


%%
function ce = kernel_ce(Prm, Blk, du, dv, A, alpha, beta)

edgewt = imbinarize(Prm,T);
% ce = A .* edgewt .* (Prm./Blk).^(alpha) .* (log(Blk./(Prm+eps))).^(beta);
end
