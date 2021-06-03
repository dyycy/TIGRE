function proj_BH = BH_ObjectRemapping(BHCalib, projlg)
%BH_OBJECTREMAPPING Summary of this function goes here
%   Detailed explanation goes here
% lgproj: log normalized projection (non-ideal BH projection)

% tranvers length grid
ulgd = BHCalib.bowtie.ulgd;
nCol = size(ulgd,2);

% Object Calibration LUT: [bowtie_sl, object_sl]
calibLUT = BHCalib.object.calibLUT;

% object_sl
object_sl = BHCalib.object.sl;

%% Remap non-ideal BH signal to linear object thickness
proj_BH = zeros(size(projlg));

% for 1D bowtie only
for jj = 1: nCol
    tl_v = ulgd(1, jj);
    % model the LUT for specific bowtie thickness
    proj_signal_LUT = interp1(BHCalib.bowtie.sl, calibLUT, tl_v, 'linear');
    proj_BH(:,jj,:) = interp1(proj_signal_LUT, object_sl, projlg(:,jj,:), 'linear');
end

end

