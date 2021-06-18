function [proj_lg, geo, angles] = VarianDataLoader(datafolder, varargin)
% VarianDataLoader: Loads Varian CBCT projection, geomtry and angles data
%   Optional parameter: Motion lag correction. Default True. 
%
% Load all dataset that are needed for reconstruction
% Tested on TrueBeam 2.0 and 2.7
% Date: 2021-04-02
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = '~/your_data_path/varian/2020-01-01_123456/';

%% GPU initialization
reset(gpuDevice(1));

%%
tag_ACDC =1
tag_DPS =1
tag_SC =1
tag_BH =1

%% Input Parser
% ACDC: acceleration & deceleration correction (default: true)
% DPS: Detector Point Spread correction (default: true)
% SC: Scatter Correction (default: true)
% BH: Beam Hardening correction (default: false, due to Bowtie, BH correction is not required at all.)
[tag_ACDC, tag_DPS, tag_SC, tag_BH] = parse_inputs(varargin{:});

%% Load geometry
[geo, ScanXML] = GeometryFromXML(datafolder);

%% Remove over-sampled projections due to acceleration and deceleration
thd = 0;
if(tag_ACDC)
    angular_interval = str2double(ScanXML.Acquisitions.Velocity.Text)...
        ./str2double(ScanXML.Acquisitions.FrameRate.Text);
    thd = angular_interval *0.9;
end

%% Load proj and angles
disp('Loading Proj: ')
tic
[proj, angles, airnorm] = ProjLoader(datafolder,thd);
toc
% Detector point scatter correction

disp('Proj DPS: ')
tic
if(tag_DPS)
    proj = DetectorPointScatterCorrection(proj, geo);
end
toc

%% Load blank scan
disp('Loading Blk: ')
tic
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);
toc
% Detector point scatter correction
disp('Blk DPS: ')
tic
if(tag_DPS)
    Blk = DetectorPointScatterCorrection(Blk, geo);
end
toc
%% Scatter Correction
tic
if(tag_SC)
    disp('Scatter correction onging: ')
    proj = ScatterCorrection(datafolder, Blk, BlkAirNorm, proj, airnorm, geo);
    disp('Scatter correction is completed.')
end
toc
%% Airnorm and Logarithmic Normalization
tic
proj_lg = LogNormal(proj, angles, airnorm, Blk, Sec, BlkAirNorm);
disp('Log Normalization is completed.')
toc
% remove anomolies
proj_lg = ZeroAnomoly(proj_lg); 

%% Beam Hardening correction is applied (kind of slow)
tic
if(tag_BH)
    [proj_lg, ~] = BHCorrection(datafolder, geo, ScanXML, proj_lg);
end
toc
%% Remove anomalies
proj_lg = ZeroAnomoly(proj_lg); 

%% mediant filtering along colume-orth
proj_lg = medfilt_col(proj_lg);

%% Gantry and Image Rotation correction: not required at all!!!
%{
[proj_lg, angles] = ImgOrient(proj_lg, angles);
%}

%% double to single
proj_lg = single(proj_lg);
angles = deg2rad(angles);

% imgFDK = FDK(proj_lg, geo, angles);
% BUG! in cropCBCT
% imcroped=cropCBCT(imgFDK, geo);

%% Audio signal 
%load train;
%sound(y,Fs)
end


function [tag_ACDC, tag_DPS, tag_SC, tag_BH] = parse_inputs(varargin)
% create input parser
p = inputParser;
% add optional parameters
addParameter(p,'acdc', true);
addParameter(p,'dps', true);
addParameter(p,'sc', true);
addParameter(p,'bh', false);

%execute
parse(p,varargin{:});
%extract
tag_ACDC=p.Results.acdc;
tag_DPS=p.Results.dps;
tag_SC=p.Results.sc;
tag_BH=p.Results.bh;
end
