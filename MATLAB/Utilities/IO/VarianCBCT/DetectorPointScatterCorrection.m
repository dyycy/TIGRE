function proj = DetectorPointScatterCorrection(proj, geo)
%% Detector Point Scatter Correction
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Date: 2021-03-26
% Author: Yi Du (yi.du@hotmail.com)

%% Empirical values from reference paper
% unit: cm-2
a0 = 3.43;
a1 = 0.000309703536035;
% unit: cm-1
a2 = 0.546566915157327;
a3 = 0.311272841141691;
% unit: cm-1
a4 = 0.002472148007134;
a5 = -12.6606856375944;

% unit: mm
offset=geo.offDetector;

% grid unit: mm
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);
% unit mm - > cm
us = us/10;
vs = vs/10;

%% Downsampling
% about 10 mm in axial direction
dus = downsample(us, 26);
% about 4 mm in axial direction
dvs = downsample(vs, 10);

%% Grid mesh
[uu,vv] = meshgrid(us,vs); %detector
[duu, dvv] = meshgrid(dus,dvs); %detector

%% Scatter convolution kernel
grid = sqrt(duu.^2 + dvv.^2);
hd = a0*(a1* exp(-a2 * grid)) + a3 * (exp( -a4 * ( grid - a5).^3 ));

%% 2D Convolution with downsampling and upsampling
for ii = 1:size(proj,3)
    page = interp2(uu, vv, proj(:,:,ii), duu, dvv);
    % gross scatter distribution
    sc = conv2(page, hd, 'same');
    % upsample the scatter distribution to the same grid level as the
    % measured intensity
    scpage = interp2(duu, dvv, sc, uu, vv,'spline');
    % primary = measure - scatter
    proj(:,:,ii) = proj(:,:,ii) - scpage;
end

%% Over correction
proj(proj<0) = eps;

end
