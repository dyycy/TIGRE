
datafolder = 'E:\BigData\RJ_CBCTscan\CatPhantom504\Inserts_2_preferred';

[geo, ScanXML] = GeometryFromXML(datafolder);

angular_interval = str2double(ScanXML.Acquisitions.Velocity.Text)...
    ./str2double(ScanXML.Acquisitions.FrameRate.Text);
thd = angular_interval *0.95;

%% load raw data
[proj, angles, airnorm] = ProjLoader(datafolder,thd);
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);

%% Data preparation
% degree to rad
theta = deg2rad(angles);

%% Part 1: reconstruction with no correction at all
raw_proj = proj;
raw_Blk = Blk;

% Log normalization
raw_proj_lg = LogNormal(raw_proj, angles, airnorm, raw_Blk, Sec, BlkAirNorm);

% in case of abnormlies
raw_proj_lg(isnan(raw_proj_lg)) = 0;
raw_proj_lg(isinf(raw_proj_lg)) = 0;

% all negative to zeros
raw_proj_lg(raw_proj_lg<0) = 0;

% mediant filtering along colume-orth
raw_proj_lg = medfilt_col(raw_proj_lg);

% Gantry Rotation correction
% Clockwise
if(angles(end) - angles(1)>0)
    raw_proj_lg = flip(raw_proj_lg, 3);
% Counter-clockwise -> Clockwise
else
    theta = flip(theta);
end

%% Part 1: recon
% raw_proj
% raw_Blk
tmp_proj = single(raw_proj_lg);
raw_img=FDK(tmp_proj,geo,theta);


%% Part 2: with detector scatter correction
% Detector point scatter correction
dsc_proj = DetectorPointScatterCorrection(proj, geo);
% Detector point scatter correction
dsc_Blk = DetectorPointScatterCorrection(Blk, geo);

% Log normalization
dsc_proj_lg = LogNormal(dsc_proj, angles, airnorm, dsc_Blk, Sec, BlkAirNorm);

% in case of abnormlies
dsc_proj_lg(isnan(dsc_proj_lg)) = 0;
dsc_proj_lg(isinf(dsc_proj_lg)) = 0;

% all negative to zeros
dsc_proj_lg(dsc_proj_lg<0) = 0;

% mediant filtering along colume-orth
dsc_proj_lg = medfilt_col(dsc_proj_lg);

% Gantry Rotation correction
% Clockwise
if(angles(end) - angles(1)>0)
    dsc_proj_lg = flip(dsc_proj_lg, 3);
% Counter-clockwise -> Clockwise
else
    theta = flip(theta);
end

%% Part 2: with detector scatter correction
tmp_proj = single(dsc_proj_lg);
dsc_img=FDK(tmp_proj,geo,theta);


%% Part 3: with object scatter correction

% Detector point scatter correction
dsc_proj = DetectorPointScatterCorrection(proj, geo);
% Detector point scatter correction
dsc_Blk = DetectorPointScatterCorrection(Blk, geo);

sc_proj = ScatterCorrection(datafolder, dsc_Blk, BlkAirNorm, dsc_proj, airnorm, geo);

sc_proj_lg = LogNormal(sc_proj, angles, airnorm, dsc_Blk, Sec, BlkAirNorm);
disp('Log Normalization is completed.')

% in case of abnormlies
sc_proj_lg(isnan(sc_proj_lg)) = 0;
sc_proj_lg(isinf(sc_proj_lg)) = 0;

% all negative to zeros
sc_proj_lg(sc_proj_lg<0) = 0;

% mediant filtering along colume-orth
sc_proj_lg = medfilt_col(sc_proj_lg);

% Gantry Rotation correction
%{ 
Clockwise
if(angles(end) - angles(1)>0)
    sc_proj_lg = flip(sc_proj_lg, 3);
% Counter-clockwise -> Clockwise
else
    theta = flip(theta);
end
%}
%% Part 3: with detector scatter correction
tmp_proj = single(sc_proj_lg);
sc_img=FDK(tmp_proj,geo,theta);


%% Part 4: with BH correction
% Key calibration information
BHCalib = BHCalibFromXML(datafolder, ScanXML);
% Precompute bowtie attenuated spectra
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
% Build bowtie attenuated spectra LUT
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
% Build reference object (water) attanuation LUT
BHCalib = BH_ObjectCalibLUT(BHCalib);


sc_proj_lg = LogNormal(sc_proj, angles, airnorm, dsc_Blk, Sec, BlkAirNorm);

% all negative to zeros
sc_proj_lg(sc_proj_lg<0) = 0;

% mediant filtering along colume-orth
sc_proj_lg = medfilt_col(sc_proj_lg);

% BH correction via reference object (water)
bh_proj_lg = BH_ObjectRemapping(BHCalib, sc_proj_lg);

% in case of abnormlies
bh_proj_lg(isnan(bh_proj_lg)) = 0;
bh_proj_lg(isinf(bh_proj_lg)) = 0;

% all negative to zeros
bh_proj_lg(bh_proj_lg<0) = 0;

% Gantry Rotation correction
% Clockwise
if(angles(end) - angles(1)>0)
    bh_proj_lg = flip(bh_proj_lg, 3);
% Counter-clockwise -> Clockwise
else
    theta = flip(theta);
end

%% Part 4: reconstruction
tmp_proj = single(bh_proj_lg);
bh_img=FDK(tmp_proj,geo,theta);


%%
slice_idx = 40;

width = 0.03;

figure(1)
subplot(2,2,1)
imshow(raw_img(:,:, slice_idx),[0 width]), title('raw')
subplot(2,2,2)
imshow(dsc_img(:,:, slice_idx),[0 width]), title('dsc')
subplot(2,2,3)
imshow(sc_img(:,:, slice_idx),[0 width]), title('sc')
subplot(2,2,4)
imshow(sc_img(:,:, slice_idx) - dsc_img(:,:,slice_idx),[-width width]/20), title('diff')


subplot(2,2,4)
imshow(bh_img(:,:, slice_idx),[0 2]), title('sc+BH')

line_idx = 256;
figure(2)
plot(raw_img(line_idx, :, slice_idx),'r'), hold on
plot(dsc_img(line_idx, :, slice_idx),'k--'), hold on
plot(sc_img(line_idx, :, slice_idx),'g'), hold on
% plot(bh_img(line_idx, :, slice_idx),'b')

%%
%{
diameter = 500;
raw_img = ROIMask(raw_img, diameter);
raw_img(raw_img<0) =0;
filename = 'raw_inserts';
write_nrrd(filename, raw_img, geo);

diameter = 500;
dsc_img = ROIMask(dsc_img, diameter);
dsc_img(dsc_img<0) =0;
filename = 'dsc_inserts';
write_nrrd(filename, dsc_img, geo);

diameter = 500;
sc_img = ROIMask(sc_img, diameter);
sc_img(sc_img<0) =0;
filename = 'sc_inserts';
write_nrrd(filename, sc_img, geo);

diameter = 500;
bh_img = ROIMask(bh_img, diameter);
bh_img(bh_img<0) =0;
filename = 'bh_inserts';
write_nrrd(filename, bh_img, geo);

%}
