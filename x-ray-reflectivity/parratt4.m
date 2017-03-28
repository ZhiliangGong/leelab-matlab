function [ Ref ] = parratt4(ed_profile, ddlay, qz, p_buff)
%Function to recursively calculate parratt reflectivity from a given
%electron density profile in units of electrons/cubic angstrom with first
%column including z coordinates and second electron density values
%Greg Tietjen Feb 17, 2012
%updated on Mar 10, 2012
%updated on Mar 15, 2012 - changed to full ed_profile input including z
%positions so that thickness can be calculated and adjusted output so that it is
%Ref/Rf
%updated on Oct 7, 2012 to calculate qc from edprofile rather than as an
%input parameter
%Daniel Kerr Mar 21, 2016 - complete overhaul


%qz = zeros(5,1);
rho = [0, ed_profile(:,2)', p_buff];
rho = rho/p_buff;
ddlay = [0, ddlay', 0];
qc = Qc_calc(p_buff);

r_fres = (sqrt(qz.^2-ones(size(qz))*rho(1)*qc^2) - sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2))./...
    (sqrt(qz.^2-ones(size(qz))*rho(1)*qc^2) + sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2));

R_fres = double(r_fres.*conj(r_fres));

R = (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) - sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2))./...
    (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) + sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2));

qj = sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2);

for j = (length(rho)-1):-1:2
    
    qjp1 = sqrt(qz.^2-ones(size(qz))*rho(j-1)*qc^2);
    reff = (qjp1 - qj)./(qjp1 + qj);
    phase = exp(1i*qj*ddlay(j));
    n1 = R.*phase;
    R = (reff + n1)./(1 + reff.*n1);
    
    if abs(rho(j-1)-rho(1)) < 1e-7
        j = 0;
    end
    
    qj = qjp1;   
    
end


Ref(:,1) = qz;
Ref(:,2) = double(R.*conj(R));
Ref(:,3) = Ref(:,2)./R_fres;

end