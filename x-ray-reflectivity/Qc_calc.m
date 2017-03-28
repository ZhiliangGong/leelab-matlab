function [ Qc ] = Qc_calc( rho )
%Function to calculate Qc given the density of the subphase
%Where rho is the desnity of the bulk and lambda the wavelength of the 
%Greg Tietjen 09/28/2011
%upated 10/02/2012 - corrected was previously not accurate

lambda = 1.24;
ro = 2.818*10^-5;% radius of electron in angstrom

k = 2*pi/lambda;

delta = 2*pi*ro*rho/k/k;

Qc = 2*k*sqrt(2*delta);



end