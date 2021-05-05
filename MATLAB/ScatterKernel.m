function [outputArg1,outputArg2] = ScatterKernel(geo,inputArg2)
%SCATTERKERNEL Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = geo;


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

%% Grid mesh
% unit: mm
[uu,vv] = meshgrid(us,vs); %detector
[duu, dvv] = meshgrid(dus,dvs); %detector


end



function ce = kernel_ce(Prm, Blk, du, dv, A, alpha, beta)

edgewt = imbinarize(Prm,T);
% ce = A .* edgewt .* (Prm./Blk).^(alpha) .* (log(Blk./(Prm+eps))).^(beta);
end
