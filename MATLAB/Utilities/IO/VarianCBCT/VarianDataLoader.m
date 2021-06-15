function [proj_lg, geo, angles] = VarianDataLoader(datafolder, varargin)
% VarianDataLoader: Loads Varian CBCT projection, geomtry and angles data
%   Optional parameter: Motion lag correction. Default True. 
%
% Load all dataset that are needed for reconstruction
% Tested on TrueBeam 2.0 and 2.7
% Date: 2021-04-02
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = '~/your_data_path/varian/2020-01-01_123456/';

%% Beam Hardening correction is applied (kind of slow)
BH = false;

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

% remove anomolies
proj_lg = ZeroAnomoly(proj_lg); 

% mediant filtering along colume-orth
% proj_lg = medfilt_col(proj_lg);

%% Beam Hardening Correction (refer to MIRT toolkit): to debug
if(BH == true)
    disp('Beam Hardening Correction is on-going: be patient... ');
    % Key calibration information
    BHCalib = BHCalibFromXML(datafolder, ScanXML);
    % Precompute bowtie attenuated spectra
    BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
    % Build reference object (water) attanuation LUT
    BHCalib = BH_ObjectCalibLUT(BHCalib);
    % BH correction via reference object (water)
    proj_lg = BH_ObjectRemapping(BHCalib, proj_lg);
    disp('BH correction is done.')
end
%% Remove anomalies
proj_lg = ZeroAnomoly(proj_lg); 

%% Gantry and Image Rotation correction
[proj_lg, angles] = ImgOrient(proj_lg, angles);

%% double to single
proj_lg = single(proj_lg);

% imgFDK = FDK(proj_lg, geo, angles);
% BUG! in cropCBCT
% imcroped=cropCBCT(imgFDK, geo);

%% Audio signal 
%load train;
%sound(y,Fs)
end
