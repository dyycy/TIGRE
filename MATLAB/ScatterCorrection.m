function prim = ScatterCorrection(datafolder, Blk, BlkAirNorm, proj, airnorm, geo)
%SCATTERKERNEL Summary of this function goes here
%   Detailed explanation goes here

% unit: mm
offset=geo.offDetector;
% Center Kernel
offset = [0, 0];

% grid unit: mm
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% Downsampling
% unit: mm
% Downsampled to about 10 mm in axial direction
dus = downsample(us, 26);
% Downsampled to about 4 mm in transaxial direction
dvs = downsample(vs, 10);

step_du = mean(diff(dus));
step_dv = mean(diff(dvs));

%% Grid
% unit: mm
[ugd,vgd] = meshgrid(us,vs); %detector
[dugd, dvgd] = meshgrid(dus,dvs); %detector

%% Load Scatter Calibration
sccalib = ScCalibFromXML(datafolder);

%% Load Blk (comment out later)
%{
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);
% Detector point scatter correction
Blk = DetectorPointScatterCorrection(Blk, geo);
%}
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

%% Scatter Correction
% Iteration Number
niter = 10;

% k(y): anti-scatter grid
ASG = ASGkernel(sccalib, geo, dus, dvs);
%% --------------------- gamma in scatter estimation
%% To debug

gamma = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{1}.ObjectScatterFit.gamma.Text);

gamma = gamma/10;

prim = zeros(size(proj));

lambda = 0.01;

for ii = 1:1 %size(proj, 3)
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
        % estimate thickness
        thickness = ThicknessEstimator(blk, page, sccalib, step_du, step_dv);
        % smooth thickness
        thickness = SmoothThickness(thickness, sccalib, step_du, step_dv);
        % Ri(x,y): group-based masks
        nmask = GroupMask(thickness, ngroup, nbounds);
        % gi(x,y): group-based form function
        gform = FormFunc(sccalib, dugd, dvgd);
        % cei(x,y): group-based amplitude factors
        cfactor = AmplitudeFactor(blk, page, sccalib);

%{
        edgewt = imbinarize(thickness, 'adaptive');
        h = fspecial('average', [1, 3]);
        edgewt = imfilter(edgewt, h);
        edgewt = imfilter(edgewt, h);
%}

        %% n-group summation
        tmp1 = 0;
        tmp2 = 0;
        for kk = 1: ngroup
            %% 2D fft
            tmp1 = tmp1 + fft2(page.*nmask(:,:,kk).*cfactor(:,:,kk)) .* fft2(gform(:,:,kk).*ASG);
            %% 2D fft
            tmp2 = tmp2 + fft2(thickness.*page.*nmask(:,:,kk).*cfactor(:,:,kk)) .* fft2(gform(:,:,kk).*ASG);            
        end
        %% -------------------------  probably complex number ------------------
        tmp1 = real(ifft2(tmp1));
        tmp2 = real(ifft2(tmp2));
        %% estimated scatter map
        Is = (1 - gamma.*thickness).*tmp1 + gamma.*tmp2;
        % Ip = Ip + lambda * (Is_prv - Is)
        page = page + lambda * (Is_prv - Is);
        page(page<0)=eps;
    end
    
    %% Upsampling and cutoff for over-correction
    % measured intensity
    SF = interp2(dugd, dvgd, Is, ugd, vgd, 'linear', 0);
    [s,l] = bounds(SF(:))
    SF(SF<0) = eps;
    SF = min(SF./proj(:,:,ii), 0.95);
    % primary
    prim(:,:,ii) = proj(:,:,ii).*(1 - SF);

end


end
