function out = Arls_nol_20110519(x,data,dflag)

% Objective function to determine Ar decay parameters
%
% x is K, lambdae, lambdab, T1...T10
% data is a structure containing the t's and r's
% dflag triggers diagnostic output

if nargin < 3; dflag = 0; end;

% Unpack
K = x(1);
le = x(2);
lb = x(3);
Ti = x(4:end)';

l = lb + le;

% Compute R's
Ri = (le./(K.*l)).*(exp(Ti.*l)-1);
% Or compute t's
% Ri = Ti;
% Ti = (1./(le+lb)).*log(((le+lb)./le).*K.*Ri + 1);

% Match to t's and r's
SS1 = 0.5*((Ri-data.ri)./data.delri).^2;
SS2 = 0.5*((Ti - data.ti)./data.delti).^2;
% Match to other constraints
%SS3 = (([K le lb l]' - data.a)./data.dela).^2;
% Ignore lambda total
SS3 = (([K le lb]' - data.a(1:3))./data.dela(1:3)).^2;



% Assemble 
out = sum(SS1)  + sum(SS2) + sum(SS3);

% Diagnostics: returns R,t,K,le,lb,l misses
if dflag == 1;
    SS1a = (Ri-data.ri)./data.delri;
    SS2a = (Ti - data.ti)./data.delti;
    SS3a = ([K le lb l]' - data.a)./data.dela;
    out = [SS1a; SS2a; SS3a];
end;