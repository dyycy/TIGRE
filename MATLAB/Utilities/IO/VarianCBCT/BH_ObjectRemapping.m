function proj_BH = BH_ObjectRemapping(BHCalib, projlg)
%BH_OBJECTREMAPPING Summary of this function goes here
%   Detailed explanation goes here
% lgproj: log normalized projection (non-ideal BH projection)

% tranvers length grid
ulgd = BHCalib.bowtie.ulgd;
[nRow, nCol] = size(ulgd);

% Object Calibration LUT: [bowtie_sl, object_sl]
calibLUT = BHCalib.object.calibLUT;

% object_sl
object_sl = BHCalib.object.sl;

%% Remap non-ideal BH signal to linear object thickness
%{
proj_BH = zeros(size(projlg));

% for 1D bowtie only
for jj = 1: nCol
    tl_v = ulgd(1, jj);
    % model the LUT for specific bowtie thickness
    proj_signal_LUT = interp1(BHCalib.bowtie.sl, calibLUT, tl_v, 'spline');
    proj_BH(:,jj,:) = interp1(proj_signal_LUT, object_sl, projlg(:,jj,:), 'spline');
end
%}

%% CPU Version: very very slow
%{
proj_BH = zeros(size(projlg));

for ii = 1: nRow
    if(~mod(ii,50))
        disp(ii);
    end
    for jj = 1: nCol
        tl_v = ulgd(ii, jj);
        % model the LUT for specific bowtie thickness
        proj_signal_LUT = interp1(BHCalib.bowtie.sl, calibLUT, tl_v, 'spline');
        proj_BH(ii,jj,:) = interp1(proj_signal_LUT, object_sl, projlg(ii,jj,:), 'spline');
    end
end
%}

%% Parfor CPU version: very slow
%{
proj_BH = zeros(size(projlg));

sl = BHCalib.bowtie.sl;
tic
for ii = 1: nRow
    if(~mod(ii,50))
        disp(ii);
        toc
        tic
    end
    parfor jj = 1: nCol
        tl_v = ulgd(ii, jj);
        % model the LUT for specific bowtie thickness
        proj_signal_LUT = interp1(sl, calibLUT, tl_v, 'spline');
        proj_BH(ii,jj,:) = interp1(proj_signal_LUT, object_sl, projlg(ii,jj,:), 'spline');
    end
end
%}

%% GPU version: acceptable
% GPU Reset
g = gpuDevice(1);
reset(g);

% gpuArray preparation
% proj_BH = single(zeros(size(projlg)));
% gproj_BH = zeros(size(proj_BH), 'gpuArray');
gprojlg = gpuArray(single(projlg));

gsl = gpuArray(BHCalib.bowtie.sl);
gLUT = gpuArray(calibLUT);

gobj_sl = gpuArray(object_sl);
proj_signal_LUT = gpuArray(zeros(1, size(gLUT,2)));

% 1D interpolation pixel-by-pixel
for ii = 1: nRow
    if(~mod(ii,50))
        disp([num2str(ii) '/' num2str(nRow)]);
    end
    for jj = 1: nCol
        tl_v = ulgd(ii, jj);
        % model the LUT for specific bowtie thickness: BUG! in interp1 for
        % spline interpolation
        % proj_signal_LUT = interp1(gsl, gLUT, tl_v, 'spline');
        proj_signal_LUT = interp1(gsl, gLUT, tl_v);
        % gproj_BH(ii, jj, :) = interp1(proj_signal_LUT, object_sl, gprojlg(:,jj,:), 'spline');
        gprojlg(ii, jj, :) = interp1(proj_signal_LUT, gobj_sl, gprojlg(ii,jj,:));
    end
end

proj_BH = double(gather(gprojlg));

% GPU reset
g = gpuDevice(1);
reset(g);

end

