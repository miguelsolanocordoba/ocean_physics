function mld = compute_mld(rho,z,crit)
% COMPUTE_MLDT Computes mixed layer depth using temperature criteria
% MLD = COMPUTE_MLDT(TEMP,CRIT) computes the mixed layer depth (MLD) from
% the Nx1 temperature profile TEMP at depths Z using tempertature criteria
% CRIT.
%
% TEMP and Z must have the same dimensions. CRIT is usually anywhere from
% 0.1C to 0.5C. Lower temperature criteria captures changes to shorter
% length/time scales (diurnal cycling) while longer temperature length
% scales capture longer length/time scales (annual cycling). 
%
% Created: May 2, 2018 by Miguel Solano

% Compute gradient of depth
gz=gradient(z); 

% Size
[nt,nz]=size(rho);

% From time 1 to nt
mld=zeros(1,nt); 
for i=1:nt
    
    grho=gradient(rho(i,:))'; 
    mflag = false;
    j=0; 
    while ~mflag
    j=j+1;
    
    % If did not find gradient > criteria
    if (j>nz)
        mflag=true; 
%         [~,ind]=max(grho./gz);
        mld(i)=mld(i-1);
    % If gradient > than criteria 
    elseif abs(grho(j)/gz(j)) > crit
        mld(i)=z(j); 
        mflag=true;
    end
    
    end
end