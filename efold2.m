function beta = efold2(x,z)
%EFOLD2 Compute e-folding depth for non-decreasing profiles. 
% 
% BETA = EFOLD2(X,Z) computes the e-folding depth BETA of variable X defined
% at vertical grid Z. X is an MxN matrix and Z is a Mx1 (or 1xM) vector, so
% that M denotes space and N denotes time. The e-folding value is the value
% at which X(z) = X(IND)*exp(-1), and [~,IND]=max(X(:,i)) at t=i. 
%
% Example: Compute e-folding depth (beta) of Stokes drift vector (us,vs); 
% 
% x = sqrt(us.^2 + vs.^2); 
% z = Vgrid; 
% beta = efold2(x,z); 
% 
% Created: May 8, 2018 by M. Solano. 

% Size of depth/time vectors
tt=size(x,2); 

 
depth=zeros(tt,1); 
for i = 1:tt
    [~,ind1]=max(x(:,i)); 
    x2=x(ind1:end,:); 
    xf = x2(1,i).*exp(-1);
    [~,ind2]=min(abs(x2(:,i)-xf));
    depth(i)=z(ind1+ind2-1);
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

%save('beta.mat') 
