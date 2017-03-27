classdef Reflectivity < handle
    
    properties
        
        rawdata
        input
        offset
        
        lipidFit
        proteinFit
        
        ed
        
    end
    
    methods
        
        function this = Reflectivity(p_buff, energy, dataFile)
            % assumes the ed profiles to be in a sub-folder called
            % 'ed-profiles'
            
            this.input.p_buff = p_buff;
            this.input.energy = energy;
            
            if nargin == 2
                [filename, pathname, ~] = uigetfile('.txt', 'Select the raw reflectivity file.');
                dataFile = fullfile(pathname, filename);
            end
            
            this.rawdata = importdata(dataFile);
            this.calculateQc();
            
        end
        
        function calculateQc(this)
            
            ro = 2.818*10^-5; % radius of electron in angstrom
            k = 2 * pi / this.wavelength();
            delta = 2 * pi * ro * this.input.p_buff / k^2;
            this.input.qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        % qz offset fitting
        
        function qzOffsetFitLipidOnly(this)
            
            qzoff_naive = this.getNaiveQzOffset();
            qc_cutoff_ind = find(this.rawdata(:, 1) <= 0.026);
            
            x0 = [qzoff_naive, 13, 10.4, 0.223, 0.453]; %qz, l_tail, l_head, p_tail, p_head
            lb = [-Inf, 0, 0, 0, 0];
            ub = [Inf, 100, 30, 100, 100];
            ProtFlag = 0;
            qzoff_edprofile = [];
            
            qzoff_fit_func = @(x, xdata) this.Qz_calc(x, xdata, qzoff_edprofile, this.input.p_buff, 3.4, this.input.qc, ProtFlag, this.rawdata);
            
            [x, resnorm] = lsqcurvefit(qzoff_fit_func, x0, this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1),...
                zeros(size(this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1))), lb, ub);
            
            qzoff_min = x(1);
            this.offset.chisqrd = resnorm;
            
            [ this.offset.ref, this.offset.refnorm, this.offset.refnormcut, ~, ~] = this.Ref_reduce3_man(this.rawdata, 0.026, qzoff_min, this.input.qc);
            
            this.input.qzoffset = qzoff_min;
            this.offset.ref = this.rawdata;
            this.offset.ref(:,1) = this.offset.ref(:,1) - this.input.qzoffset;
            
            fprintf('\n%s %f \n','The Qz offset is:', this.input.qzoffset)
            
        end
        
        function qzOffsetFitLipidProtein(this)
            
            this.importEdProfiles();
            
            qzoff_naive = this.getNaiveQzOffset();
            qc_cutoff_ind = find(this.rawdata(:, 1) <= 0.026);
            
            x0 = [qzoff_naive, -10, 0.35, 13, 10.4, 0.223, 0.453]; %qz, prot_pos, cov, l_tail, l_head, p_tail, p_head
            lb = [-Inf, -100, 0, 0, 0, 0, 0];
            ub = [Inf, 0, 1, 100, 30, 100, 100];
            ProtFlag = 1;
            
            qzoff_edprofile = this.ed.profiles{1, 1};
            
            qzoff_fit_func = @(x, xdata) this.Qz_calc(x, xdata, qzoff_edprofile, this.input.p_buff, 3.4, this.input.qc, ProtFlag, this.rawdata);
            
            [x, resnorm] = lsqcurvefit(qzoff_fit_func, x0, this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1),...
                zeros(size(this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1))), lb, ub);
            
            qzoff_min = x(1);
            this.offset.chisqrd = resnorm;
            
            [ this.offset.ref, this.offset.refnorm, this.offset.refnormcut,~,~] = ...
                this.Ref_reduce3_man(this.rawdata,0.026,qzoff_min,this.input.qc);
            
            this.input.qzoffset = qzoff_min;
            this.offset.ref = this.rawdata;
            this.offset.ref(:,1) = this.offset.ref(:,1) - this.input.qzoffset;
            fprintf('\n%s %f \n','The Qzoffset is:', this.input.qzoffset)
            
        end
        
        % data fitting
        
        function dataFitLipid(this, sigma)
            
            if nargin == 1
                sigma = 3.4;
            end
            
            this.qzOffsetFitLipidOnly();
            
            fprintf('\n %s\n', 'The current version cuts all data below 0.026 for fitting purposes - this data is typically unreliable');
            cutind = find(this.offset.refnorm(:,1) > 0.026, 1, 'first');
            this.lipidFit.refnorm = this.offset.refnorm(cutind+1:end,:);
            
            x0 = [13, 10.4, 0.223, 0.453]; %l_tail, l_head, p_tail, p_head
            lb = [0, 0, 0, 0];
            ub = [100, 30, 100, 100];
            
            lipid2box_func = @(x,xdata) this.lipid2box(x, xdata, this.lipidFit.refnorm(:,3), this.input.p_buff, 0, sigma);
            [x, this.lipidFit.chi2] = lsqcurvefit(lipid2box_func, x0, this.lipidFit.refnorm(:,1), this.lipidFit.refnorm(:,2)./this.lipidFit.refnorm(:,3), lb, ub);
            
            this.lipidFit.l_tail = x(1);
            this.lipidFit.l_head = x(2);
            this.lipidFit.p_tail = x(3);
            this.lipidFit.p_head = x(4);
            this.lipidFit.sigma = sigma;
            
            this.lipidFit.p_buff = this.input.p_buff;
            
            [this.lipidFit.ed, ddlay] = this.Lipid_Prot_EDcalc4([], [], [], x(1), x(2), x(3), x(4), this.lipidFit.p_buff, this.lipidFit.sigma, 0);
            this.lipidFit.ref_fit = this.parratt4(this.lipidFit.ed, ddlay, this.lipidFit.refnorm(:,1), this.lipidFit.p_buff);
            
        end
        
        function dataFitProtein(this, sigma)
            
            if nargin == 1
                sigma = 3.4;
            end
            
            this.qzOffsetFitLipidProtein();
            
            fprintf('\n %s\n', 'The current version cuts all data below 0.026 for fitting purposes - this data is typically unreliable');
            cutind = find(this.offset.refnorm(:, 1) > 0.026, 1, 'first');
            this.proteinFit.refnorm = this.offset.refnorm(cutind+1:end, :);
            
            x0 = [-10, 0.35, 13, 10.4, 0.223, 0.453]; %prot_pos, cov, l_tail, l_head, p_tail, p_head
            lb = [-100, 0, 0, 0, 0, 0];
            ub = [0, 1, 100, 30, 100, 100];
            
            theta = this.ed.theta;
            phi = this.ed.phi;
            edProfiles = this.ed.profiles;
            refnorm = this.proteinFit.refnorm;
            m = length(theta);
            n = length(phi);
            d = length(x0);
            x = zeros(m, n, d);
            chi = zeros(m, n);
            p_buff = this.input.p_buff;
            opts = optimset('Display', 'off');
            
            tic;
            parfor i = 1 : m
                refnorm_par = refnorm;
                for j = 1 : n
                    prot_ed = edProfiles{i, j};
                    prot2box_func = @(x, xdata) prot2box(x, xdata, refnorm_par(:, 3), prot_ed, p_buff, 0, sigma);
                    [p, chi_value] = lsqcurvefit(prot2box_func, x0, refnorm_par(:,1), refnorm_par(:,2) ./ refnorm_par(:,3), lb, ub, opts);
                    x(i, j, :) = reshape(p, 1, 1, length(p));
                    chi(i, j) = chi_value;
                end
            end
            toc;
            
            this.proteinFit.para = x;
            this.proteinFit.chi = chi;
            
            this.proteinFit.position = x(:, :, 1);
            this.proteinFit.coverage = x(:, :, 2);
            this.proteinFit.tailLength = x(:, :, 3);
            this.proteinFit.headLength = x(:, :, 4);
            this.proteinFit.tailEd = x(:, :, 5);
            this.proteinFit.headEd = x(:, :, 6);
            this.proteinFit.bufferEd = p_buff;
            this.proteinFit.sigma = sigma;
            
        end
        
        % utility
        
        function importEdProfiles(this)
            
            pattern = '^t\d\d\dp\d\d\d.ed$';
            path = uigetdir(pwd, 'Select the folder for ED files.');
            
            try
                contents = dir(path);
                if isempty(contents)
                    path = uigetdir(pwd, 'Select the folder for ED files.');
                    this.importEdProfiles();
                end
            catch
                this.importEdProfiles();
            end
            
            isEdFile = false(1, length(contents));
            allFiles = cell(1, length(contents));
            for i = 1 : length(contents)
                allFiles{i} = contents(i).name;
                if regexp(contents(i).name, pattern) == 1
                    isEdFile(i) = true;
                end
            end
            
            this.ed.path = path;
            this.ed.files = allFiles(isEdFile);
            
            names = repmat('0', length(this.ed.files), length(this.ed.files{1}));
            for n = 1 : length(this.ed.files)
                names(n, :) = this.ed.files{n};
            end
            
            theta = sort(str2num(names(:, 2 : 4)));
            phi = sort(str2num(names(:, 6 : 8)));
            
            theta = theta(theta ~= [theta(end); theta(1 : end -1)]);
            phi = phi(phi ~= [phi(end); phi(1 : end -1)]);
            
            if length(theta) * length(phi) ~= length(this.ed.files)
                error('Some ed files are missing, or their names are incorrect');
            else
                this.ed.theta = theta;
                this.ed.phi = phi;
                
                this.ed.profiles = cell(length(theta), length(phi));
                disp('Importing all the ed files, will take a moment...');
                for n = 1 : length(theta)
                    for m = 1 : length(phi)
                        filename = ['t', sprintf('%03d', theta(n)), 'p', sprintf('%03d', phi(m)), '.ed'];
                        file = fullfile(path, filename);
                        fid = fopen(file);
                        this.ed.profiles{n, m} = cell2mat(textscan(fid, '%f %f %f'));
                        fclose(fid);
                    end
                end
            end
            
        end
        
        function qzoff = getNaiveQzOffset(this)
            
            % Determine data points below qc:
            
            ind = this.rawdata(:, 2) > 0.95 * this.rawdata(1, 2);
            norm = mean(this.rawdata(ind, 2));
            
            qz = this.rawdata(:, 1);
            refdata = this.rawdata(:,2) / norm;
            
            % Define fit function
            % ft = fittype('Rfcalc3(x,qc,qzoff)','problem','qc');
            Rfcalcfit = @(qzoff, xdata) this.calculateReflectivity(xdata, this.input.qc, qzoff);
            
            %Perform fit and extract qzoff
            %f = fit(qz,refdata,ft, 'problem',qc,'StartPoint', .0004);
            %qzoff = coeffvalues(f);
            qzoff = lsqcurvefit(Rfcalcfit, 0.0004, qz, refdata, [], []);
            
        end
        
        function wavelength = wavelength(this)
            %wavelength in A
            
            c = 299792458;
            h = 6.62607004e-34;
            ev2j = 1.60218e-19;
            
            wavelength = h * c / (this.input.energy * 1000 * ev2j) * 1e10;
            
        end
        
        function R = Qz_calc(this, x, reflxaxis, prot_ed, p_buff, sigma, qc, ProtFlag, refl)
            % Function for fitting qz offset to 2 box model with/without protein
            
            [~, ~, Refl_Fresnel_qz_shift_cutoff, ~, ~] = this.Ref_reduce3_man(refl, 0.026, x(1), qc);
            
            if ProtFlag == 1
                [ED, ddlay] = this.Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, sigma, ProtFlag);
            elseif ProtFlag == 0
                [ED, ddlay] = this.Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), p_buff, sigma, ProtFlag);
            end
            
            Reflparratt = this.parratt4(ED, ddlay, reflxaxis-x(1), p_buff);
            
            R = (Refl_Fresnel_qz_shift_cutoff(:,2) - Reflparratt(:,3))./Refl_Fresnel_qz_shift_cutoff(:,3);
            
            
        end
        
        function Ref = parratt4(this, ed_profile, ddlay, qz, p_buff)
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
            rho = rho / p_buff;
            ddlay = [0, ddlay', 0];
            qc = this.input.qc;
            
            r_fres = (sqrt(qz.^2 - ones(size(qz)) * rho(1) * qc^2) - sqrt(qz.^2 - ones(size(qz)) * rho(end) *qc^2)) ./ ...
                (sqrt(qz.^2 - ones(size(qz)) * rho(1) * qc^2) + sqrt(qz.^2 - ones(size(qz)) * rho(end) * qc^2));
            
            R_fres = double(r_fres .* conj(r_fres));
            
            R = (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) - sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2))./...
                (sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2) + sqrt(qz.^2-ones(size(qz))*rho(end)*qc^2));
            
            qj = sqrt(qz.^2-ones(size(qz))*rho(end-1)*qc^2);
            
            for j = (length(rho)-1) : -1 : 2
                
                qjp1 = sqrt(qz.^2-ones(size(qz))*rho(j-1)*qc^2);
                reff = (qjp1 - qj)./(qjp1 + qj);
                phase = exp(1i*qj*ddlay(j));
                n1 = R .* phase;
                R = (reff + n1)./(1 + reff.*n1);
                
                qj = qjp1;
                
            end
            
            
            Ref(:,1) = qz;
            Ref(:,2) = double(R.*conj(R));
            Ref(:,3) = Ref(:,2)./R_fres;
            
        end
        
        function [R] = lipid2box(this, x, reflxaxis, reflerr, p_buff, sigmaflag, sigma)
            
            %Function for fitting lipid 2box model to qz offset reflectivity data
            
            if sigmaflag == 1
                [ED, ddlay] = this.Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), p_buff, x(1) , 0);
            elseif sigmaflag == 0
                [ED, ddlay] = this.Lipid_Prot_EDcalc4([], [], [], x(1), x(2), x(3), x(4), p_buff, sigma, 0);
            end
            
            Reflparratt = this.parratt4(ED, ddlay, reflxaxis, p_buff);
            
            R = Reflparratt(:,3)./reflerr;
            
        end
        
        function [R] = prot2box(this, x, reflxaxis, reflerr, prot_ed, p_buff, sigmaflag, sigma)
            
            %Function for fitting lipid 2box model to qz offset reflectivity data
            
            if sigmaflag == 1
                [ED, ddlay] = this.Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, x(1) , 1);
            elseif sigmaflag == 0
                [ED, ddlay] = this.Lipid_Prot_EDcalc4(prot_ed, x(1), x(2), x(3), x(4), x(5), x(6), p_buff, sigma, 1);
            end
            
            Reflparratt = this.parratt4(ED, ddlay, reflxaxis, p_buff);
            
            R = Reflparratt(:,3)./reflerr;
            
        end
        
        % plot
        
        function plotLipidFit(this)
            
            figure;
            errorbar(this.lipidFit.refnorm(:, 1), this.lipidFit.refnorm(:, 2), this.lipidFit.refnorm(:, 3), '.', 'color', [0, 0, 0] + 0.5, 'linewidth', 1.2);
            hold on;
            plot(this.lipidFit.ref_fit(:, 1), this.lipidFit.ref_fit(:, 3), '-k', 'linewidth', 2.4);
            set(gca, 'fontsize', 14);
            hold off;
            xlabel('$$ Q_z (\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('Normalized Reflectivity', 'fontsize', 16);
            legend('Data', 'Fit');
            
        end
        
    end
    
    methods(Static)
        
        function Rf = calculateReflectivity(qz_obs, qc, qzoff)
            
            qz = qz_obs-qzoff;

            r_fres = (qz - sqrt(qz.^2 - ones(size(qz)) * qc^2)) ./ (qz + sqrt(qz.^2 - ones(size(qz))*qc^2));
            Rf = double(r_fres.*conj(r_fres));
            
        end
        
        function [ oRef, oRefoRf, ocRefoRf, qc, qzoff] = Ref_reduce3_man(Ref, qc_cut, qzoff, qc)
            %Function for processing Raw Reflectivity Data
            
            format long;
            
            %Determine data points below qc:
            
            ind = find(Ref(:,1)<=qc_cut);
            
            
            %Rf generation
            
            
            Rf = Rfcalc3(Ref(:,1), qc, qzoff);
            
            oRef = Ref;
            oRef(:,1) = oRef(:,1)-qzoff;
            oRefoRf(:,1) = Ref(:,1)-qzoff;
            oRefoRf(:,2) = Ref(:,2)./Rf;
            oRefoRf(:,3) = Ref(:,3)./Rf;
            
            
            %Rf cut generation
            
            ocRefoRf = oRefoRf([2:4 (length(ind)+1):end],:);
            
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
            
        end
        
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
        
    end
    
end