function [h,Rib] = compute_bld(u,v,T,S,z,t,ustar,B0,f,Ric)
%COMPUTE_BLD computes BLD based on bulk Richardson number (Large et al. 1994).
% 
% [H,RIB] = compute_bld(u,v,T,S,z,t,ustar,B0,f,Ric) computes boundary layer 
% depth  H based on the bulk Richardson number RIB, as computed in 
% Large et al. (1994). 
% 
% u = velocity in the x-direction (nz,nt) [m/s]
% v = velocity in the y-direction (nz,nt) [m/s]
% T = temperature (nz,nt) [C]
% S = salinity (nz,nt) [PSU]
% z = depths (nz) [m]
% t = time (nt) 
% ustar = friction velocity (nt) [m/s]
% B0 = surface buoyancy flux (nt) [m^2/s^3]
% f = Coriolis frequency [1/s]
% Ric = Critical Richardson number
%
% The bulk Richardson number is the ratio of buoyancy to a velocity 
% referenced from the surface. 
%
% Rib = (B_r - B_r(z))*z / (|V_r - V_r(z)| + V_t(z))
%
% Refer to Large et al. (1994) and Smyth et al. (2002) for details. 
%
% Created: Feb 19, 2019 by M. Solano
% Modified: Aug 7, 2019 by M. Solano

%% Constants
g = 9.81;      % Gravity 
rho0 = 1025;   % Seawater density (constant)
kappa = 0.4;   % Von Karman constant
betaT = -0.2;  % ratio of entrainment flux to surface buoyancy flux
epsilon = 0.1; % Surface layer fraction 
Cv = 1.6;      % calibrated to yield beta=-0.2 under pure convection
as = -28.86;   % coefficient of nondimensional flux profile for scalar
cs = 98.96;    % coefficient of nondimensional flux profile for scalar
phis = -1.0;   % Non-dimensional flux

nz = numel(z); % number of vertical layers
nt = numel(t); % number of time steps

d = abs(z);    % Distance from surface (always positive) 

%% T/S derived quantities 
% Compute potential density and stratification frequency
pden = sw_pden(S,T,d,0);
[~,dRhodz] = gradient(pden,t,z);
N2 = -g*dRhodz./pden;
% N2 = -g*dRhodz./rho0;

% Velocity and buoyancy profiles
V = sqrt(u.^2 + v.^2); 
B = -g.*pden./rho0; 

% Mixed layer depth, Ekman depth and Obukhov length
%hm = compute_mld(T,S,z,0.125); 
hm = compute_mld2(T',z,0.1); 
he = -0.7*ustar/f; 
L = ustar.^3./(kappa.*B0); 

%% Compute turbulent velocity for scalars (ws) 
% Non-dimensional stability parameter (zeta) and vertical coordinate
% (sigma). Note that the vertical coordinate is normalized using the MLD
% instead of BLD. 
zeta = d./L; 
sigma = d./abs(hm);

% Nondimensional flux(phis) [for details see Large et al. 1994, Appendix B]
phi_s = (1 - 16*zeta).^(-1/2);
phi_s(zeta>=0) = 1 + 5.*zeta(zeta>=0);
phi_s(zeta<phis) = (as - cs*zeta(zeta<phis)).^(-1/3);
ws = kappa*ustar./phi_s; 

% For unstable conditions (zeta<0) ws(sigma>0.1)=ws(sigma=0.1)
for i = 1:nt
    for j = 1:nz
        if zeta(j,i)<0 && sigma(j,i)>0.1
            [~,ind] = min(abs(sigma(:,i)-0.1)); 
            ws(j,i) = kappa*ustar(i)/phi_s(ind,i);
        else
            ws(j,i) = kappa*ustar(i)/phi_s(j,i); 
        end
    end
end


% turbulent velocity scale [Large et al. 1994, equation (23)]
Vt2 = Cv*(-betaT)^0.5/(Ric*kappa^2) * (cs*epsilon)^(-0.5) * (d.*real(sqrt(N2)).*ws);

% pre-allocation
Rib = zeros(nz,nt); 
h = zeros(1,nt); 
ht = zeros(1,nt+1); 
ht(1) = hm(1); 


%% Estimate BLD. 
[~,inds] = min(abs(10+z)); 
for i = 1:nt
    [~,ind] = min(abs(ht(i)*epsilon - z));
    ind(ind<2)=2; 
    Br = trapz(abs(z(1:ind)),B(1:ind,i))/abs(z(ind)-z(1));
    Vr = trapz(abs(z(1:ind)),V(1:ind,i))/abs(z(ind)-z(1));
    zsl = (z(ind)+z(1))/2;
    Rib(:,i) = ((Br-B(:,i)).*(zsl-z))./((Vr-V(:,i)).^2+Vt2(:,i)); 
    indt = max([inds ind]);
    [~,indh] = min(abs(Rib(indt:end,i) - Ric)); 
    h(i) = z(indt+indh-1);
    ht(i+1) = h(i); 
    if (L(i)>0 && h(i)<he(i))
        h(i) = he(i); 
        ht(i+1) = h(i); 
    end
end

save('Ribprofile.mat','Rib')
