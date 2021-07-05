function BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib)
%SPECTRAPROCESS Summary of this function goes here
%   Detailed explanation goes here

%% Thin window filter has been applied to the spectrum already
% Source-to-detector distance: mm
SDD = geo.DSD;
% Source-to-bowtie distance: mm
SBD = BHCalib.bowtie.distance;

%% [u,v] vector: mm
offset=geo.offDetector;
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% to bowtier surface
ub = us * SBD/SDD;
% vb = vs * SBD/SDD;

u_s = fix(min(ub)) + sign(min(us));
u_l = fix(max(ub)) + sign(max(us));

sampling_rate = 1002;

xxgd = linspace(u_s, u_l, sampling_rate);
y_l = ceil(max(BHCalib.bowtie.thickness));
yygd = linspace(-y_l, y_l, sampling_rate);
zzgd = linspace(-20, 20, 502);

bt_img = zeros(length(xxgd), length(yygd));

%% Fitting: first fitting(linear)
[xData, yData] = prepareCurveData(BHCalib.bowtie.uu, BHCalib.bowtie.thickness);
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[~, gof] = fit( xData, yData, ft );
if( gof.rsquare < 0.995)
    disp('It seems that BH correction is not required at all');
    typein = input('Continue BH correction? Y or N (recommended): ', 's');
    if(contains(typein, 'n', 'IgnoreCase', true))
        BHCalib = NaN;
        return;
    end
end

%% Fitting: Smoothing Spline 
% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft);
disp(['Bowtie Fitting Goodness (R^2): ' num2str(gof.rsquare)])

%% 3D bowtie img matrix
for ii = 1:length(xxgd)
    tmp_thick = fitresult(xxgd(ii));
    tmp = ((yygd>0)&(yygd<tmp_thick));
    bt_img(tmp, ii) = 1;
end

bt_img = flipud(bt_img);

bt_img3D = single(repmat(bt_img, [1 1 length(zzgd)]));

%% geo for bowtie projection
geo_bt = geo;
geo_bt.DSO = BHCalib.bowtie.distance;
geo_bt.dVoxel = [mean(diff(yygd)); mean(diff(xxgd)); mean(diff(zzgd))];
geo_bt.nVoxel = size(bt_img3D)';
geo_bt.sVoxel = geo_bt.dVoxel.* geo_bt.nVoxel;

%% tranverse length grid: mm
ulgd = Ax(bt_img3D, geo_bt, 0, 'interpolated');

%%
%{
% tranverse length vector: mm
ul = interp1(BHCalib.bowtie.uu, BHCalib.bowtie.thickness, ub);
% tranverse length grid: mm
ulgd = repmat(ul, [length(vs), 1]);
%}

%% bowtie material attenuation look-up table
[Min, Max] = bounds(ulgd(:));
% sampling length
sl = linspace(Min, Max, 101);
% LUT: [sl, miu_bowtie]
atten_tab = zeros(length(sl), length(BHCalib.bowtie.ac));

% I/I_0 = exp(- length * miu): dependent on incident energy
for ii = 1:length(sl)
    atten_tab(ii,:) = exp(-sl(ii).* BHCalib.bowtie.ac');
end

% kVp
kVp = length(BHCalib.filter.spec);
% miu
BHCalib.object.ac =BHCalib.object.ac(1:kVp);
%% bowtie-attenuated spectra look-up table
% specLUT(sampling length, energy bin)
BHCalib.bowtie.specLUT = atten_tab(:,1:kVp).*BHCalib.filter.spec;

%% sampling length
BHCalib.bowtie.sl = sl;

%% tranvers length of each detector unit
BHCalib.bowtie.ulgd = ulgd;

end
