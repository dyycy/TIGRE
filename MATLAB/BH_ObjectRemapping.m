function [outputArg1,outputArg2] = BH_ObjectRemapping(BHCalib, lgproj)
%BH_OBJECTREMAPPING Summary of this function goes here
%   Detailed explanation goes here
% lgproj: log normalized projection (non-ideal BH projection)

% 
ulgd = BHCalib.bowtie.ulgd;
[nRow, nCol] = size(ulgd);
ul_vs = reshape(ulgd, [nRow*nCol 1]);

% Object Calibration LUT
calibLUT = BHCalib.object.calibLUT;
% LUT grid for 2D interpolation
% [x, y]: object_sl, bowtie_sl
xgd = BHCalib.object.xgdLUT;
ygd = BHCalib.object.ygdLUT;

% 
for ii = 1:nRow
    for jj = 1: nCol
        bowtie_thickness = ulgd(ii,jj);
        object_thickness = 

        proj_signal = lgproj(ii,jj);

        
%% just shift y, x
%%%%%%%%%%%%%%%        
for mm=1:MM
	s1hat(:,mm) = interp1(f1(:,mm), s1, fh1(:,mm), 'pchip', 'extrap');
end
%%%%%%%%%%%%%%%

end

