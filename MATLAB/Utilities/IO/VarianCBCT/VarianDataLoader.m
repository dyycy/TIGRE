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
disp('Scatter correction onging: ')
proj_sc = ScatterCorrection(datafolder, Blk, BlkAirNorm, proj, airnorm, geo);
disp('Scatter correction is completed.')

%% Airnorm and Logarithmic Normalization
proj_lg = LogNormal(proj_sc, angles, airnorm, Blk, Sec, BlkAirNorm);
disp('Log Normalization is completed.')

% all negative to zeros
proj_lg(proj_lg<0) = 0;

% mediant filtering along colume-orth
proj_lg = medfilt_col(proj_lg);

%% Beam Hardening Correction: refer to MIRT toolkit
% Key calibration information
BHCalib = BHCalibFromXML(datafolder, ScanXML);
% Precompute bowtie attenuated spectra
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
% Build bowtie attenuated spectra LUT
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
% Build reference object (water) attanuation LUT
BHCalib = BH_ObjectCalibLUT(BHCalib);
% BH correction via reference object (water)
proj_BH = BH_ObjectRemapping(BHCalib, projlg);

%% Mediate filtering along colume-orth
proj = medfilt_col(proj);

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
img = FDK(proj,geo,angles);

end


%% Mediant filtering along colume-orth
%{
function proj = medfilt_col(proj)

% 3D projection matrix
for ii = 1:size(proj,3)
    proj(:,:,ii) = ordfilt2(proj(:,:,ii), 5, ones(1,9));
end

end
%}

