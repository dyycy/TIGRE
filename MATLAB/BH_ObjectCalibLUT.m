function BHCalib = BH_ObjectCalibLUT(BHCalib)
%BHWATERCALIB Summary of this function goes here
%   Detailed explanation goes here

% max object thickness: mm
MaxThickness = BHCalib.object.thicknessmax;

% object sampling length: mm
object_sl = linspace(0, MaxThickness, round(MaxThickness));

%% spectra[sI, energy]
specLUT = BHCalib.bowtie.specLUT;
miu = BHCalib.object.ac;

% object thickness look-up table: [bowtie.sl, object_sl]
calibLUT = zeros(length(BHCalib.bowtie.sl), length(object_sl));

% bowtie thickness - > spectrum
for ii = 1: length(BHCalib.bowtie.sl)    
    % object thickness -> 
    spec = specLUT(ii,:);
    for jj = 1:length(object_sl)
        % non-ideal attenuated signal
        tmp = spec .* exp(-object_sl(jj).* miu);
        % nonlinear BH projection
        calibLUT(ii,jj) = -log(sum(tmp)/sum(spec));
    end
end

% non ideal projection signal LUT: [bowtie.sl, object_sl]
BHCalib.object.calibLUT = calibLUT;

% grids for interpolation: bowtie thickness, obgd
% be careful with the X-Y orientation in meshgrid
[obj_gd, bt_gd] = meshgrid(object_sl, BHCalib.bowtie.sl);

% grid for LUT
BHCalib.object.xgdLUT = obj_gd;
BHCalib.object.ygdLUT = bt_gd;

end

