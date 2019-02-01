function beta = efold(x,z)
%EFOLD Compute e-folding depth. 
% 
% BETA = EFOLD(X,Z) computes the e-folding depth BETA of variable X defined
% at vertical grid Z. X is an MxN matrix and Z is a Mx1 (or 1xM) vector, so
% that M denotes space and N denotes time. The e-folding value is the value
% at which X(z) = X(1)*exp(-1). X is assumed to be maximum at X(1),
% otherwise use EFOLD2!
%
% Example: Compute e-folding depth (beta) of Stokes drift vector (us,vs); 
% 
% x = sqrt(us.^2 + vs.^2); 
% z = Vgrid; 
% beta = efold(x,z); 
% 
% Created: May 8, 2018 by M. Solano. 

% Size of depth/time vectors
tt=size(x,2); 

xf = x(1,:).*exp(-1); 
depth=zeros(tt,1); 
for i = 1:tt
    [~,ind]=min(abs(x(:,i)-xf(i)));
    depth(i)=z(ind);
%     if x(ind,i)-xf(i) > 0
%         deltaz=z(ind+1)-z(ind); 
%         deltax=x(ind+1,i)-x(ind,i); 
%         deltaxn=xf(i)-x(ind,i); 
%         depth(i)=z(ind)+deltaxn*deltaz/deltax;
%     elseif x(ind,i)-xf(i) < 0 
%         deltaz=z(ind)-z(ind-1); 
%         deltax=x(ind,i)-x(ind-1,i);
%         deltaxn=xf(i)-x(ind-1,i);
%         depth(i)=z(ind-1)+deltaxn*deltaz/deltax; 
%     else
%         depth(i) = z(ind);
%     end
end

beta=1./depth; 

% save('beta.mat') 
