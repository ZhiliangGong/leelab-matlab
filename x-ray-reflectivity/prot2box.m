function [R] = prot2box(x, reflxaxis, reflerr, prot_ed, p_buff, sigmaflag, sigma)

%Function for fitting lipid 2box model to qz offset reflectivity data

if sigmaflag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, x(1) , 1);
elseif sigmaflag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(1), x(2), x(3), x(4), x(5), x(6), p_buff, sigma, 1);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis, p_buff);

R = Reflparratt(:,3) ./ reflerr;

end