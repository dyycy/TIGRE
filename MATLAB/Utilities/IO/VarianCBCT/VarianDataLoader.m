function [proj,geo, angles] = VarianDataLoader(datafolder, varargin)
% VarianDataLoader: Loads Varian CBCT projection, geomtry and angles data
%   Optional parameter: Motion lag correction. Default True. 
%
% Load all dataset that are needed for reconstruction
% Tested on TrueBeam 2.0 and 2.7
% Date: 2021-04-02
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = '~/your_data_path/varian/2020-01-01_123456/';

%% Load geometry
[geo, ScanXML] = GeometryFromXML(datafolder);

%% Remove the over-sampled projections due to acceleration and decceleration
thd = 0;
if(~isempty(varargin)&&(varargin{1}))
    angular_interval = str2double(ScanXML.Acquisitions.Velocity.Text)...
        ./str2double(ScanXML.Acquisitions.FrameRate.Text);
    thd = angular_interval *0.95;
end

%% Load proj and angles
[proj, angles, airnorm] = ProjLoader(datafolder,thd);
% Detector point scatter correction
proj = DetectorPointScatterCorrection(proj, geo);

%% Load blank scan
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);
% Detector point scatter correction
Blk = DetectorPointScatterCorrection(Blk, geo);

%% Scatter Correction
proj = ScatterCorrection(datafolder, Blk, BlkAirNorm, proj, airnorm, geo);

%% Airnorm and Logarithmic Normalization
proj = LogNormal(proj, angles, airnorm, Blk, Sec, BlkAirNorm);

%% Beam Hardening Correction via MIRT toolkit



%% Mediate filtering along colume-orth
for ii = 1:size(proj,3)
    proj(:,:,ii) = ordfilt2(proj(:,:,ii), 5, ones(1,9));
end

% in case of abnormlies
proj(isnan(proj)) = 0;
proj(isinf(proj)) = 0;

% all negative to zeros
proj(proj<0) = 0;

% double to single
proj = single(proj);

% degree to rad
angles = angles/180*pi;

%% Gantry Rotation correction
% Clockwise
if(angles(end) - angles(1)>0)
    proj = flip(proj, 3);
% Counter-clockwise -> Clockwise
else
    angles = flip(angles);
end


%%

% img=FDK(proj,geo,angles);

end
