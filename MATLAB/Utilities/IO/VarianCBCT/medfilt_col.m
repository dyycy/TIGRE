function proj = medfilt_col(proj)


% 3D projection matrix
for ii = 1:size(proj,3)
    proj(:,:,ii) = ordfilt2(proj(:,:,ii), 5, ones(1,9));
end


end


