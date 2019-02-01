function mld = compute_mldt2(temp,z,crit)
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

% Compute gradients
gtemp=gradient(temp)'; 
gz=gradient(z); 

mflag = false;
i=0; 
numel(gtemp);
max(abs(gtemp./gz));
while ~mflag
    i=i+1
    if (i>numel(gtemp)) && (max(abs(gtemp./gz))>0.1)
        [~,ind]=max(gtemp./gz); 
        mflag=true; 
        mld=z(ind); 
    elseif (i>numel(gtemp)) && (max(abs(gtemp./gz))<=0.1)
        mflag=true; 
        mld=z(end); 
    elseif abs(gtemp(i)/gz(i)) > crit
        mld=z(i); 
        mflag=true;
    end
end