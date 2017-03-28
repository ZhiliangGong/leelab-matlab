function Rf = calculateReflectivity(qz_obs, qc, qzoff)

qz = qz_obs - qzoff;

r_fres = (qz - sqrt(qz.^2 - ones(size(qz)) * qc^2)) ./ (qz + sqrt(qz.^2 - ones(size(qz))*qc^2));
Rf = double(r_fres.*conj(r_fres));

end