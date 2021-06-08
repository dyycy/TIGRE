function proj = ZeroAnomoly(proj)
%% Remove anomalies
% in case of NaN or Inf
proj(isnan(proj)) = 0;
proj(isinf(proj)) = 0;

% all negative to zeros
proj(proj<0) = 0;


end

