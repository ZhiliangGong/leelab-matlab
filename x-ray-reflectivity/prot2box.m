function [R] = prot2box(x, reflxaxis, reflerr, prot_ed, p_buff, sigmaflag, sigma)

%Function for fitting lipid 2box model to qz offset reflectivity data

if sigmaflag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, x(1) , 1);
elseif sigmaflag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(1), x(2), x(3), x(4), x(5), x(6), p_buff, sigma, 1);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis, p_buff);

R = Reflparratt(:,3) ./ reflerr;

    function [ ED_out, ddlay ] = Lipid_Prot_EDcalc4( prot_ed, protpos, cov, l_tail, l_head, p_tail, p_head, p_buff, sigma, ProtFlag)
        %Greg Tietjen Oct 7 2012 - Most current version for use
        %Daniel Kerr Mar 7 2016 - Added contingency if protein does not penetrate
        %tail group at all
        %Daniel Kerr Mar 15 2015 - Complete overhaul
        
        %Program to generate Lipid/Protein electron density profile in z dimension
        %given parameters derived from x ray scattering fit with Schlossman method.
        
        
        L_lipid = l_head + l_tail;
        
        %First need to properly orient prot_ed position values. They are initially
        %in a different coordinate system - first value => protpos
        
        if ProtFlag == 1
            prot_z = prot_ed(:,1) - prot_ed(1,1) + protpos;
            thickness = prot_z(end-1) - prot_z(end);
            %     prot_z = zeros(size(prot_z0));
            %
            %     for i = 1:length(prot_z)
            %         prot_z(i) = prot_z0(1) - thickness*(i-1);
            %     end
            
            prot_ed_vac = prot_ed(:,2);
            prot_emp_vol = prot_ed(:,3);
            
            %Next determine how much of protein is in lipid layer vs how much is in
            %buffer
            
            z_profile = [0, prot_z'];
            %z_profile0 = [0, prot_z0'];
            %thickness = abs(prot_z(1,1)-prot_z(2,1));
            %dlay = thickness*ones(size(prot_z(:,1)));
            %dlay = [protpos, dlay];
            e_dens(1) = p_tail;
            
            %Define interface indices from z_profile
            
            tail_head_ind = find(z_profile >= -l_tail,1,'last');
            tailheadindepro = find(prot_z >= -l_tail,1,'last');
            
            %Check if protein extends into tail region
            
            if isempty(tailheadindepro) == 0
                %Make sure that bottom of protein is not at the tail-head interface
                if tail_head_ind < length(z_profile)
                    %Check if l_tail is already listed in protein ed profile, if not,
                    %add index to ed profile to include interface position and assign
                    %it the same density as the layer before, but correct the thickness
                    %of that layer
                    if z_profile(tail_head_ind)~= -l_tail
                        
                        z_profile = [z_profile(1:tail_head_ind), -l_tail, z_profile(tail_head_ind+1:end)];
                        
                        tail_head_ind = tail_head_ind + 1;
                        %Shift index to now refer to actual l_tail position
                        
                        %dlay = [dlay(1:tail_head_ind-1), (z_profile(tail_head_ind-1)-z_profile(tail_head_ind)),...
                        %    (z_profile(tail_head_ind)-z_profile(tail_head_ind+1)), dlay(tail_head_ind, end)];
                        e_dens(2:tail_head_ind-1) = cov*prot_ed_vac(1:tailheadindepro) + (1+cov*(prot_emp_vol(1:tailheadindepro)-1))*p_tail;
                        e_dens(tail_head_ind) = cov*prot_ed_vac(tailheadindepro) + (1+cov*(prot_emp_vol(tailheadindepro)-1))*p_head;
                        
                    else
                        e_dens(2:tail_head_ind) = cov*prot_ed_vac(1:tailheadindepro) + (1+cov*(prot_emp_vol(tailheadindepro)-1))*p_tail;
                    end
                    
                    prof_start = tail_head_ind + 1;
                    prot_start = tailheadindepro + 1;
                    
                else
                    
                    if z_profile(tail_head_ind)~= -l_tail
                        
                        z_profile = [z_profile(1:tail_head_ind), -l_tail];
                        
                        tail_head_ind = tail_head_ind + 1;
                        %Shift index to now refer to actual l_tail position
                        
                        %dlay = [dlay(1:tail_head_ind-1), (z_profile(tail_head_ind-1)-z_profile(tail_head_ind))];
                        e_dens(2:(length(prot_ed_vac)+1)) = cov*prot_ed_vac(1:end) + (1+cov*(prot_emp_vol(1:end)-1))*p_tail;
                        e_dens((length(prot_ed_vac)+1)) = p_tail;
                        
                    else
                        e_dens(2:(length(prot_ed_vac)+1)) = cov*prot_ed_vac(1:end) + p_tail;
                    end
                    
                    %dlay = [dlay, l_head];
                    e_dens = [e_dens, p_head];
                    prof_start = tail_head_ind + 1;
                    prot_start = tailheadindepro + 1;
                    
                end
                
            else
                z_profile = [0, -l_tail, z_profile(2:end)];
                e_dens(2) = p_head;
                prof_start = 3;
                prot_start = 1;
            end
            
            
            head_buffer_ind = find(z_profile >= -L_lipid,1,'last');
            headbufferindepro = find(prot_z >= -L_lipid,1,'last');
            
            if isempty(headbufferindepro) == 0
                %Check bottom of protein is not above head-buffer interface
                if head_buffer_ind < length(z_profile)
                    %Check if L_lipid is already listed in protein ed profile, if not,
                    %add index to ed profile to include interface position and assign
                    %it the same density as the layer before, but correct the thickness
                    %of that layer
                    if z_profile(head_buffer_ind)~= -L_lipid
                        
                        z_profile = [z_profile(1:head_buffer_ind), -L_lipid, z_profile(head_buffer_ind+1:end)];
                        
                        head_buffer_ind = head_buffer_ind+1;
                        %Shift index to now refer to actual L_lipid position
                        
                        e_dens(prof_start:head_buffer_ind-1) = cov*prot_ed_vac(prot_start:headbufferindepro) + (1+cov*(prot_emp_vol(prot_start:headbufferindepro)-1))*p_head;
                        e_dens(head_buffer_ind) = cov*prot_ed_vac(headbufferindepro) + (1+cov*(prot_emp_vol(headbufferindepro)-1))*p_buff;
                        
                    else
                        e_dens(prof_start:head_buffer_ind) = cov*prot_ed_vac(prot_start:headbufferindepro) + (1+cov*(prot_emp_vol(prot_start:headbufferindepro)-1))*p_head;
                    end
                    
                else
                    
                    if z_profile(head_buffer_ind)~= -L_lipid
                        
                        z_profile = [z_profile(1:head_buffer_ind), -L_lipid];
                        
                        head_buffer_ind = head_buffer_ind+1;
                        %Shift index to now refer to actual L_lipid position
                        
                        e_dens(prof_start:head_buffer_ind-1) = cov*prot_ed_vac(prot_start:headbufferindepro) + (1+cov*(prot_emp_vol(prot_start:headbufferindepro)-1))*p_head;
                        e_dens(head_buffer_ind) = p_head;
                        
                    else
                        e_dens(prof_start:head_buffer_ind) = cov*prot_ed_vac(prot_start:end) + (1+cov*(prot_emp_vol(prot_start:end)-1))*p_head;
                    end
                    
                end
                
                prof_start = head_buffer_ind + 1;
                prot_start = headbufferindepro + 1;
                
                if head_buffer_ind < length(z_profile)
                    e_dens(prof_start:(length(prot_ed_vac)-(prot_start))+prof_start) =...
                        cov*prot_ed_vac(prot_start:end) + (1+cov*(prot_emp_vol(prot_start:end)-1))*p_buff;
                end
                
            elseif isempty(headbufferindepro) == 1 && isempty(tailheadindepro) == 1
                
                z_profile = [z_profile(1:2), -L_lipid, z_profile(3:end)];
                e_dens = [p_tail, p_head, p_buff, (cov*prot_ed_vac' + (1+cov*(prot_emp_vol'-1))*p_buff)];
                
                %prof_start = 4;
                %prot_start = 1;
                
            end
            
            
            z_profile = [z_profile, (z_profile(end)-thickness)];
            e_dens = [0, e_dens, p_buff];
            
        elseif ProtFlag == 0
            
            z_profile = [0, -l_tail, -L_lipid];
            e_dens = [0, p_tail, p_head, p_buff];
            
        end
        
        
        %Next need to smooth using error function at each interface
        
        if ProtFlag == 1
            
            if sigma > 0
                
                ed_length = (z_profile(1) - z_profile(end) + 14*sigma)/0.5;
                ed_step = (z_profile(1) - z_profile(end) + 14*sigma)/(ed_length+1);
                
                ddlay = zeros(floor(ed_length+1), 1);
                ED_out = zeros(floor(ed_length+1), 2);
                
                for ed_i = 1:(ed_length+1)
                    ddlay(ed_i) = ed_step;
                    za =  z_profile(1) + 7*sigma - ed_i*ed_step;
                    chem = 0;
                    %z_min = z_profile(1) - za;
                    %chem = chem + (0/2)*erfc(z_min/(sqrt(2)*sigma));
                    
                    for j = 1:length(z_profile)
                        
                        z_min = z_profile(j) - za;
                        
                        if j == 1
                            chem = chem + (e_dens(j)/2)*erf(z_min/(sqrt(2)*sigma));
                        else
                            z_max = z_profile(j-1) - za;
                            chem = chem + (e_dens(j)/2)*(erf(z_max/(sqrt(2)*sigma))-erf(z_min/(sqrt(2)*sigma)));
                        end
                        
                    end
                    
                    z_max = z_profile(end) - za;
                    chem = chem + (e_dens(length(z_profile)+1)/2)*(1+erf(z_max/(sqrt(2)*sigma)));
                    %chem = chem/p_buff ;
                    ED_out(ed_i,1) = za;
                    ED_out(ed_i,2) = chem;
                    
                end
                
            else
                
                ed_length = (z_profile(1) - z_profile(end) + 20)/0.5;
                ed_step = (z_profile(1) - z_profile(end) + 20)/(ed_length+1);
                
                ddlay = zeros(floor(ed_length+1), 1);
                ED_out = zeros(floor(ed_length+1), 2);
                
                for ed_i = 1:(ed_length+1)
                    ddlay(ed_i) = ed_step;
                    za =  z_profile(1) + 10 - ed_i*ed_step;
                    chem = 0;
                    %z_min = z_profile(1) - za;
                    %chem = chem + (0/2)*erfc(z_min/(sqrt(2)*sigma));
                    
                    if za > z_profile(end) && za < z_profile(1)
                        for j = 1:length(z_profile)
                            if za >= z_profile(j) && za < z_profile(j-1)
                                chem=e_dens(j);
                            end
                        end
                    elseif za >= z_profile(1)
                        chem = 0;
                    else
                        chem = e_dens(length(z_profile)+1);
                    end
                    
                    %chem = chem/p_buff ;
                    ED_out(ed_i,1) = za;
                    ED_out(ed_i,2) = chem;
                    
                end
                
            end
            
        elseif ProtFlag == 0
            
            if sigma > 0
                
                ed_length = (z_profile(1) - z_profile(end) + 14*sigma)/0.25;
                ed_step = (z_profile(1) - z_profile(end) + 14*sigma)/(ed_length+1);
                
                ddlay = zeros(floor(ed_length+1), 1);
                ED_out = zeros(floor(ed_length+1), 2);
                
                for ed_i = 1:(ed_length+1)
                    ddlay(ed_i) = ed_step;
                    za =  z_profile(1) + 7*sigma - ed_i*ed_step;
                    chem = 0;
                    
                    for j = 1:length(z_profile)
                        
                        z_min = z_profile(j) - za;
                        if j == 1
                            chem = chem + (e_dens(j)/2)*erf(z_min/(sqrt(2)*sigma));
                        else
                            z_max = z_profile(j-1) - za;
                            chem = chem + (e_dens(j)/2)*(erf(z_max/(sqrt(2)*sigma))-erf(z_min/(sqrt(2)*sigma)));
                        end
                        
                    end
                    
                    z_max = z_profile(end) - za;
                    chem = chem + (e_dens(length(z_profile)+1)/2)*(1+erf(z_max/(sqrt(2)*sigma)));
                    %chem = chem/p_buff ;
                    ED_out(ed_i,1) = za;
                    ED_out(ed_i,2) = chem;
                    
                end
                
            else
                
                ed_length = (z_profile(1) - z_profile(end) + 20)/0.25;
                ed_step = (z_profile(1) - z_profile(end) + 20)/(ed_length+1);
                
                ddlay = zeros(floor(ed_length+1), 1);
                ED_out = zeros(floor(ed_length+1), 2);
                
                for ed_i = 1:(ed_length+1)
                    ddlay(ed_i) = ed_step;
                    za =  z_profile(1) + 10 - ed_i*ed_step;
                    chem = 0;
                    %z_min = z_profile(1) - za;
                    %chem = chem + (0/2)*erfc(z_min/(sqrt(2)*sigma));
                    
                    if za > z_profile(end) && za < z_profile(1)
                        for j = 1:length(z_profile)
                            if za >= z_profile(j) && za < z_profile(j-1)
                                chem=e_dens(j);
                            end
                        end
                    elseif za >= z_profile(1)
                        chem = 0;
                    else
                        chem = e_dens(length(z_profile));
                    end
                    
                    %chem = chem/p_buff ;
                    ED_out(ed_i,1) = za;
                    ED_out(ed_i,2) = chem;
                    
                end
                
            end
            
        end
        
        
    end

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
        
        R_temp = (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) - sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2))./...
            (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) + sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2));
        
        qj = sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2);
        
        for j = (length(rho)-1):-1:2
            
            qjp1 = sqrt(qz.^2-ones(size(qz))*rho(j-1)*qc^2);
            reff = (qjp1 - qj)./(qjp1 + qj);
            phase = exp(1i*qj*ddlay(j));
            n1 = R_temp.*phase;
            R_temp = (reff + n1)./(1 + reff.*n1);
            
            qj = qjp1;
            
        end
        
        
        Ref(:,1) = qz;
        Ref(:,2) = double(R_temp.*conj(R_temp));
        Ref(:,3) = Ref(:,2)./R_fres;
        
    end

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

end