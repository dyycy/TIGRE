function BHCalib = BH_ObjectCalibLUT(BHCalib)
%BHWATERCALIB Summary of this function goes here
%   Detailed explanation goes here

% max object thickness: mm
MaxThickness = BHCalib.object.thicknessmax;

% object sampling length: mm
object_sl = linspace(0, MaxThickness, ceil(MaxThickness)*10);

%% spectra[sI, energy]
specLUT = BHCalib.bowtie.specLUT;
miu = BHCalib.object.ac;

% object thickness look-up table: [bowtie.sl, object_sl]
calibLUT = zeros(length(BHCalib.bowtie.sl), length(object_sl));

%% Photo flux based calibration
%{
% bowtie thickness - > spectrum
for ii = 1: length(BHCalib.bowtie.sl)    
    % object thickness -> [min_bowtie_thickness: max_bowtie_thickness]
    spec = specLUT(ii,:);
    for jj = 1:length(object_sl)
        % non-ideal attenuated signal
        tmp = spec .* exp(-object_sl(jj).* miu);
        % nonlinear BH projection
        calibLUT(ii,jj) = -log(sum(tmp)/sum(spec));
    end
end
%}

%% Energy Flux based calibration
% scintillator energy absorption coefficient
scintillator_ac = BHCalib.scintillator.ac(1: length(miu));
scintillator_thick = BHCalib.scintillator.thickness;

% bowtie thickness - > spectrum
for ii = 1: length(BHCalib.bowtie.sl)    
    % object thickness -> [min_bowtie_thickness: max_bowtie_thickness]
    spec = specLUT(ii,:);
    
    % incident energy fluence
    fluence_0 = sum(BHCalib.source.kV .*spec .*scintillator_ac.*scintillator_thick);
    
    for jj = 1:length(object_sl)
        % attenuated spectram 
        tmp = spec .* exp(-object_sl(jj).* miu);
        
        % attenuated energy fluence
        fluence_1 = sum(tmp .*BHCalib.source.kV .*scintillator_ac.*scintillator_thick);
        
        % nonlinear BH projection
        calibLUT(ii,jj) = -log(fluence_1/fluence_0);
    end
end
%% non ideal projection signal LUT: [bowtie.sl, object_sl]
BHCalib.object.calibLUT = calibLUT;

% object sampling length
BHCalib.object.sl = object_sl;

end

