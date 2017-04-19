function [ oRef, oRefoRf, ocRefoRf, qc, qzoff] = Ref_reduce3_man(Ref, qc_cut, qzoff, qc)
%Function for processing Raw Reflectivity Data

%Determine data points below qc:

ind = find(Ref(:,1)<=qc_cut);


%Rf generation


Rf = Rfcalc3(Ref(:,1),qc,qzoff);

oRef = Ref;
oRef(:,1) = oRef(:,1)-qzoff;
oRefoRf(:,1) = Ref(:,1)-qzoff;
oRefoRf(:,2) = Ref(:,2)./Rf;
oRefoRf(:,3) = Ref(:,3)./Rf;


%Rf cut generation

ocRefoRf = oRefoRf([2:4 (length(ind)+1):end],:);

end
