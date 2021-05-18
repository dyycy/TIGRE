function thickness = EstimateThickness(Blk, BlkAirNorm, Prm, AirNorm, sccalib, step_du, step_dv)
%% Estimate Water-Equivalent Thickness 
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               Blk: Total intensity, i.e., I_0
%               BlkAirNorm: 
%               Prm: Primary intensity, i.e., I_p
%               AirNorm: Primary intensity airnorm chamber value
%               sccalib: Scatter Calibration Structure
% Output:
%               thickness: Estimated object thickness, i.e., tau(x,y) in Reference
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

% mu H2O = 0.02 /mm
muH2O = str2double(sccalib.CalibrationResults.Globals.muH2O.Text);

% Air Chamber Normalization Factor
NF = AirNorm/BlkAirNorm; 

% unit mm
thickness = log(NF*Blk./Prm) /muH2O;

%% Smooth the estimated thickness
% thickness(vv, uu, ntheta)
thickness = SmoothThickness(thickness, sccalib, step_du, step_dv);

end

