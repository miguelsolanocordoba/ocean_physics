function [mld,ind] = compute_mld2(x,z,crit)
% COMPUTE_MLD2 Computes mixed layer depth based on difference from surface.
% MLD = COMPUTE_MLD2(x,z,CRIT) computes the mixed layer depth (MLD) by
% calculating the difference crit in the surface value of x and depth
% z=MLD. x Can be temperature or density, and is a 2D variable with size
% (nt,nz), where nt is the number of time samples and nz is the number of
% vertical layers. z is a one dimensional depth (nz,1) and crit is a
% scalar. For temperature crit is usually 0.2 and for density it is usually
% 0.125-0.4. 
%
% Created: August, 2018 by Miguel Solano

% Size
[nt,~]=size(x);

% Increase resolution 
z2=interp1(z,linspace(1,numel(z),500));

rho2=zeros(nt,500); 
for i=1:nt
   rho2(i,:)=interp1(z',x(i,:),z2); 
end

% From time 1 to nt
mld=zeros(1,nt);
ind=zeros(1,nt); 
for i=1:nt
   [~,ind(i)]=min(abs(rho2(i,:)-rho2(i,1)-crit)); 
   mld(i)=z2(ind(i));
end

