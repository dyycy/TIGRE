function thickness = EstimateThickness(Proj, Prm, sccalib)
%% Estimate Water-Equivalent Thickness 
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               Proj: Total intensity, i.e., I_0
%               Prm: Primary intensity, i.e., I_p
%               sccalib: Scatter Calibration Structure
% Output:
%               thickness: Estimated object thickness, i.e., tau(x,y) in Reference
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

% mu H2O = 0.02 /mm
muH2O = sccalib.CalibrationResults.Globals.muH2O;
% unit mm
thickness = log(Proj./Prm) /muH2O;

%% Remember to smooth the estimated thickness


end

