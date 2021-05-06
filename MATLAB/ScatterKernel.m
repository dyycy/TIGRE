function [outputArg1,outputArg2] = ScatterKernel(datafolder, geo,ScanXML)
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




end



function ce = kernel_ce(Prm, Blk, du, dv, A, alpha, beta)

edgewt = imbinarize(Prm,T);
% ce = A .* edgewt .* (Prm./Blk).^(alpha) .* (log(Blk./(Prm+eps))).^(beta);
end
