function BHCalibXML = BHCalibFromXML(datafolder)
%% Load Calibration.xml for Scatter Correction
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               datafolder: Varian CBCT scan data folder
% Output:
%               ScCalibXML: Scatter Calibration Structure
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

srcfilename = [datafolder filesep 'Calibrations' filesep 'ASC' filesep 'Factory' filesep 'Calibration.xml'];

%% Export as struct
tmp = xml2struct(srcfilename);
BHCalibXML = tmp.Calibration;

%% Restructure BHCalibXML: spectrum data
for ii =1:length(BHCalibXML.CalibrationResults.Spectra.SpectrumProperties)
    tmp = cell2mat(BHCalibXML.CalibrationResults.Spectra.SpectrumProperties{ii}.Flux.float);
    for jj = 1:length(tmp)
        flux(jj) = str2double(tmp(jj).Text);
    end
    % Sepctrum
    BHCalibXML.CalibrationResults.Spectra.SpectrumProperties{ii}.Spec = flux;
end

%% Restructure BHCalibXML: material attenuation coefficients


end
