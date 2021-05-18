function cfactor = AmplitudeFactor(Prj, Prm, sccalib)
%% AMPLITUDEFACTOR Summary of this function goes here

%% group number
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);

% Amplitude factor groups
cfactor = [];

for ii=1:ngroup
    tmp = sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.ObjectScatterFit;
    % Amplitude Factor
    % unit: cm^2 - > mm^2
    A = str2double(tmp.A.Text) * 100;
    % unitless
    alpha = str2double(tmp.alpha.Text);
    beta = str2double(tmp.beta.Text);
    cfactor(:,:,ii) = A.* (Prm./(Prj+eps)).^(alpha) .* ( log(Prj./(Prm+eps)) ).^(beta);    
end

end

