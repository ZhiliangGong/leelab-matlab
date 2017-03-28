function [Rf] = Rfcalc3(qz_obs,qc,qzoff)
%Rfcalc - program to calculate Fresnel reflectivity given qz and
%qc (critical edge q) and qz offset (experimentally determined error in qz)
%Greg Tietjen 04/13/2011
%Daniel Kerr 03/30/2016 Updated to compute Fresnel reflectivity based on
%electron density only


qz = qz_obs-qzoff;

r_fres = (qz - sqrt(qz.^2-ones(size(qz))*qc^2))./...
    (qz + sqrt(qz.^2-ones(size(qz))*qc^2)); 

Rf = double(r_fres.*conj(r_fres));

end