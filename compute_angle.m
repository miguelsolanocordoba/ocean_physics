function theta = compute_angle(ux,uy,vx,vy)
%COMPUTE_ANGLE compute angle between vectors (ux,uy) and (vx,vy).
% 
% THETA = COMPUTE_ANGLE(ux,uy,vx,vy) computes the angle THETA between
% vector U=(ux,uy) and V=(vx,vy). ux/ux and vx/vy are the components of
% vectors U and V respectively, in the x and y directions. Each vector
% component (ux,uy,vx,vy) must be the same size. 
%
% Example: Compute the angle between vectors Us=(us,vs) and tau=(taux,tauy). 
%
% theta = compute_angle(us(1,:),vs(1,:),taux(1,1:ind),tauy(1,1:ind)); 
% 
% Created: January 28, 2019 by M. Solano. 

% Compute magnitude of u=(ux,uy) and v = (vx,vy); 
magu = sqrt(ux.^2 + uy.^2); 
magv = sqrt(vx.^2 + vy.^2);

% Compute angle between u and v.
theta = acosd((ux.*vx + uy.*vy)./(magu.*magv));
