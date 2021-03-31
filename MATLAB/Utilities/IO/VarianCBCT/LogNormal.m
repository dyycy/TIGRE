function proj = LogNormal(proj, angles, airnorm, Blk, Sec, BlkAirNorm)
% Log Normalization: 
% Calculate Logrithmic Projections with AirNorm
% Tested on TrueBeam 2.5 and 2.7
% Date: 2021-03-23
% Author: Yi Du (yi.du@hotmail.com)

%% TB version
% Version = 2.5
if(isempty(Sec))
    for ii = 1:length(angles)
        CF = airnorm(ii)/BlkAirNorm;
        proj(:,:,ii) = log(CF*Blk./proj(:,:,ii));
    end
% Version = 2.7    
else
    % interpolation weights
    for ii = 1:length(angles)
        [left, weights] = interp_weight(angles(ii), Sec);
        interp_blk = weights(1) * Blk(:,:,left) + weights(2) * Blk(:,:,left+1);
        % Correction factor
        CF = airnorm(ii)/(0.5*BlkAirNorm(left) + 0.5*BlkAirNorm(left+1));
        proj(:,:,ii) = log(CF*interp_blk./proj(:,:,ii));
    end
end

end
