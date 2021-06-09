function [proj, angles] = ImgOrient(proj, angles)
% image reorienation
% angles output should be in radian
% Date: 2021-06-08
% Author: Yi Du, 
% Email: yi.du@hotmail.com

%% Part 1: Gantry Rotation Alignment
% Clockwise
if(angles(end) - angles(1)>0)
    proj = flip(proj, 3);
% Counter-clockwise -> Clockwise
else
    angles = flip(angles);
end

%% Part 2: image reoritation
% flag for unit
tag = (range(angles) <3*pi);

if(tag)
    disp('angles are in radian: ');
    angles = rad2deg(angles);
end

% full scan: image reorientation is not needed
if(range(angles) > 270)
    disp('Full arc scan:')
    angles = deg2rad(angles);
    return;
end


% Offset for image matrix reorientation
GantryOffset = 20;

tmp = (angles(1) + GantryOffset)/90;
% Realign to orthoganal angular position
if( abs(tmp - round(tmp)) <0.1 )
   angles = angles + GantryOffset;
end

% angles output should be in radian
angles = deg2rad(angles);

end
