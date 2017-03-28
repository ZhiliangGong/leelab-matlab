function [R] = lipid2box(x, reflxaxis, reflerr, p_buff, sigmaflag, sigma)

%Function for fitting lipid 2box model to qz offset reflectivity data

if sigmaflag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), p_buff, x(1) , 0);
elseif sigmaflag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(1), x(2), x(3), x(4), p_buff, sigma, 0);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis, p_buff);

R = Reflparratt(:,3)./reflerr;

end