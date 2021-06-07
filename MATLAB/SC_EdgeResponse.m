function edgewt = SC_EdgeResponse(thickness)
% EDGERESPONSE Summary of this function goes here
% Date: 2021-05-24

% Emperical Value
edgewt = double(imbinarize(thickness, 50));
tmpmask = edgewt;
h = fspecial('average', [25, 25]);
for ii = 1:5
    edgewt = imfilter(edgewt, h);
end        
tmp = tmpmask.*edgewt;
tmp(tmp==0) = NaN;
tmp = rescale(tmp, 0.6, 1);
tmp(isnan(tmp)) = 0;
edgewt = tmp;

end

