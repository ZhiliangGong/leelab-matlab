function R = Qz_calc(x, reflxaxis, prot_ed, p_buff, sigma, qc, ProtFlag, refl)
% Function for fitting qz offset to 2 box model with/without protein

[~, ~, Refl_Fresnel_qz_shift_cutoff, ~, ~] = Ref_reduce3_man(refl, 0.026, x(1), qc);

if ProtFlag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, sigma, ProtFlag);
elseif ProtFlag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), p_buff, sigma, ProtFlag);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis-x(1), p_buff);

R = (Refl_Fresnel_qz_shift_cutoff(:,2) - Reflparratt(:,3))./Refl_Fresnel_qz_shift_cutoff(:,3);


end