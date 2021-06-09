function prim = ScatterCorrection(datafolder, Blk, BlkAirNorm, proj, airnorm, geo)
% Scatter Correction Module
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% No downsampling or upsampling is involved in current version
% Author: Yi Du (yi.du@hotmail.com)
% Date: 2021-05-24

%% Center Coordinates
%unit: mm
offset = geo.offDetector;
offset = [0, 0];

% grid unit: mm
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% Downsampling (not used in current version)
% downsampling rate
ds_rate = 8;
% unit: mm
% Downsampled to about 10 mm in axial direction in Ref
dus = decimate(us, ds_rate);
% Downsampled to about 4 mm in transaxial direction in Ref
dvs = decimate(vs, ds_rate);

step_du = mean(diff(dus));
step_dv = mean(diff(dvs));

%% X, Y grid
% unit: mm
[ugd,vgd] = meshgrid(us,vs); %detector
% downsampled grid
[dugd, dvgd] = meshgrid(dus,dvs); %detector

%% Load Scatter Calibration
sccalib = ScCalibFromXML(datafolder);

%% Blk scan
sBlk = sum(Blk, 3);
sAirNorm = sum(BlkAirNorm);

%% n-thickness group number and boundaries
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);
nbounds = [];

for ii=1:ngroup
    % unit: mm
    tmp = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.Thickness.Text);
    nbounds = [nbounds, tmp];
end

%% k(y): anti-scatter grid kernel
ASG = SC_ASGkernel(sccalib, geo, dus, dvs);

%% Component Weights: gamma (gamma = 0 for SKS)
gamma = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{1}.ObjectScatterFit.gamma.Text);
% unit: cm-> mm
mm2cm =  1/10;
% gamma = gamma * unit_cvt;

%% iteration number
niter = 8;
% relaxation factor in iteration
lambda = 0.005;

%% Primary signal matrix
prim = zeros(size(proj));

for ii = 1: size(proj, 3)
    if(~mod(ii,50))
        display([num2str(ii) '/' num2str(size(proj, 3))]);
    end
    
    % blk: I_0 unattenuated blk signal
    CF = sAirNorm/airnorm(ii);
    blk = interp2(ugd, vgd, sBlk/CF, dugd, dvgd, 'linear', 0);       
    
    % page: I_p primary signal
    % note: detector point spread deconvolution should be done first
    page = interp2(ugd, vgd, proj(:,:,ii), dugd, dvgd, 'linear', 0);

    % Is initialization
    Is = zeros(size(page));
    
    %% Iterative Correction
    for jj = 1: niter
        % Is previous scatter map
        Is_prv = Is;
        
        % estimate thickness: mm
        thickness = SC_ThicknessEstimator(blk, page, sccalib, step_du, step_dv);
        % smooth thickness
        thickness = SC_SmoothThickness(thickness, sccalib, step_du, step_dv);
        
        % Ri(x,y): group-based masks
        nmask = SC_GroupMask(thickness, ngroup, nbounds);

        % gi(x,y): group-based form function
        gform = SC_FormFunc(sccalib, dugd, dvgd);
        
        % edge response function: about 7.5 cm inward
        edgewt = SC_EdgeResponse(thickness);
        
        % cei(x,y): group-based amplitude factors
        cfactor = SC_AmplitudeFactor(blk, page, edgewt, sccalib);
        
        % mm -> cm
        thickness = thickness * mm2cm;
        %% n-group summation
        comp1 = zeros(size(page));
        comp2 = zeros(size(page));
        for kk = 1: ngroup
            %% 2D fft
            comp1 = comp1 + fft2(page.*nmask(:,:,kk).*cfactor(:,:,kk)) .* fft2(gform(:,:,kk).*ASG);
            %% 2D fft
            comp2 = comp2 + fft2(thickness.*page.*nmask(:,:,kk).*cfactor(:,:,kk)) .* fft2(gform(:,:,kk).*ASG);            
        end
        %% real components cutoff
        comp1 = real(ifft2(comp1));
        comp2 = real(ifft2(comp2));
        %% fASKS scatter correction
        Is = (1 - gamma .*thickness).*comp1 + gamma.*comp2; 
        %sum(sum(Is_prv - Is))
        page = page + lambda * (Is_prv - Is);
        page(page<0) = eps;
    end
    
    %% Upsampling and cutoff for over-correction
    % measured intensity
    scmap = interp2(dugd, dvgd, Is, ugd, vgd, 'spline');
    % extrolation errors
    scmap(isnan(scmap)) = eps;
    % SF = Is;
    scmap(scmap<0) = eps;

    % scatter fraction
    SF = scmap./proj(:,:,ii);
    % in case of zeros in proj
    SF(SF==Inf) = NaN;
    SF(SF>1000) = NaN;
    SF = inpaint_nans(SF,0);
    % mediant filtering
    SF = medfilt2(SF, [3 3]);

    % Scatter fraction cutoff
    SF = min(SF, 0.95);    
    % primary 
    prim(:,:,ii) = proj(:,:,ii).*(1 - SF);
end

end
