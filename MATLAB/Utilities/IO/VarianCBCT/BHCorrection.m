function [proj_lg, BHCalib] = BHCorrection(datafolder, geo, ScanXML, proj_lg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('Beam Hardening Correction is on-going: be patient... ');

% Key calibration information
BHCalib = BHCalibFromXML(datafolder, ScanXML);

%->
% Precompute filter attenuated spectrum
BHCalib = BH_SpectrumFilter(geo, BHCalib);
%->

% Precompute bowtie attenuated spectra
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);
% Build reference object (water) attanuation LUT
BHCalib = BH_ObjectCalibLUT(BHCalib);
% BH correction via reference object (water)
proj_lg = BH_ObjectRemapping(BHCalib, proj_lg);

disp('BH correction is done.')

end

