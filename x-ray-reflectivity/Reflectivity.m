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
        
        % qz offset fitting
        
        function qzOffsetFitLipidOnly(this)
            
            qzoff_naive = this.getNaiveQzOffset();
            qc_cutoff_ind = find(this.rawdata(:, 1) <= 0.026);
            
            x0 = [qzoff_naive, 13, 10.4, 0.223, 0.453]; %qz, l_tail, l_head, p_tail, p_head
            lb = [-Inf, 0, 0, 0, 0];
            ub = [Inf, 100, 30, 100, 100];
            ProtFlag = 0;
            qzoff_edprofile = [];
            
            qzoff_fit_func = @(x, xdata) Qz_calc(x, xdata, qzoff_edprofile, this.input.p_buff, 3.4, this.input.qc, ProtFlag, this.rawdata);
            
            [x, resnorm] = lsqcurvefit(qzoff_fit_func, x0, this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1),...
                zeros(size(this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1))), lb, ub);
            
            qzoff_min = x(1);
            this.offset.chisqrd = resnorm;
            
            [ this.offset.ref, this.offset.refnorm, this.offset.refnormcut, ~, ~] = Ref_reduce3_man(this.rawdata, 0.026, qzoff_min, this.input.qc);
            
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
            
            qzoff_fit_func = @(x, xdata) Qz_calc(x, xdata, qzoff_edprofile, this.input.p_buff, 3.4, this.input.qc, ProtFlag, this.rawdata);
            
            [x, resnorm] = lsqcurvefit(qzoff_fit_func, x0, this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1),...
                zeros(size(this.rawdata([2:4 (length(qc_cutoff_ind)+1):end],1))), lb, ub);
            
            qzoff_min = x(1);
            this.offset.chisqrd = resnorm;
            
            [ this.offset.ref, this.offset.refnorm, this.offset.refnormcut,~,~] = ...
                Ref_reduce3_man(this.rawdata,0.026,qzoff_min,this.input.qc);
            
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
            
            lipid2box_func = @(x,xdata) lipid2box(x, xdata, this.lipidFit.refnorm(:,3), this.input.p_buff, 0, sigma);
            [x, this.lipidFit.chi2] = lsqcurvefit(lipid2box_func, x0, this.lipidFit.refnorm(:,1), this.lipidFit.refnorm(:,2) ./ this.lipidFit.refnorm(:,3), lb, ub);
            
            this.lipidFit.l_tail = x(1);
            this.lipidFit.l_head = x(2);
            this.lipidFit.p_tail = x(3);
            this.lipidFit.p_head = x(4);
            this.lipidFit.sigma = sigma;
            
            this.lipidFit.p_buff = this.input.p_buff;
            
            [this.lipidFit.ed, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(1), x(2), x(3), x(4), this.lipidFit.p_buff, this.lipidFit.sigma, 0);
            this.lipidFit.ref_fit = parratt4(this.lipidFit.ed, ddlay, this.lipidFit.refnorm(:,1), this.lipidFit.p_buff);
            
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
            
            ed_par = cell(m, n);
            ref_par = cell(m, n);
            
            tic;
            parfor i = 1 : m
                para_par = x(i, :, :);
                refnorm_par = refnorm;
                for j = 1 : n
                    prot_ed = edProfiles{i, j};
                    [ed_par{i, j}, ddlay] = Lipid_Prot_EDcalc4(prot_ed, para_par(1, j, 1), para_par(1, j, 2), para_par(1, j, 3), para_par(1 ,j, 4), para_par(1, j, 5), para_par(1, j, 6), p_buff, sigma, 1);
                    ref_par{i, j} = parratt4(ed_par{i ,j}, ddlay, refnorm_par(:, 1), p_buff);
                end
            end
            toc;
            
            this.proteinFit.ed = ed_par;
            this.proteinFit.ref_fit = ref_par;
            
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
            Rfcalcfit = @(qzoff, xdata) calculateReflectivity(xdata, this.input.qc, qzoff);
            
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
        
        function calculateQc(this)
            
            ro = 2.818*10^-5; % radius of electron in angstrom
            k = 2 * pi / this.wavelength();
            delta = 2 * pi * ro * this.input.p_buff / k^2;
            this.input.qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        % plot
        
        function plotLipidFit(this)
            
            this.plotDataAndFit(this.lipidFit.refnorm, this.lipidFit.ref_fit);
            
        end
        
        function plotChiMap(this)
            
            figure;
            chi = this.proteinFit.chi;
            theta = this.ed.theta;
            phi = this.ed.phi;
            
            imagesc(theta, phi, chi' - min(chi(:)));
            
            ylabel('\phi (deg)', 'fontsize', 16);
            xlabel('\theta (deg)', 'fontsize', 16); 
            set(gca,'YDir','normal', 'fontsize', 14);
            colorbar;
            colormap(this.pmap);
            caxis([0, chi2inv(0.999999,2)]);
            colorbar('YTick',[chi2inv(0.95,2), chi2inv(0.99,2), chi2inv(0.999,2), chi2inv(0.9999,2), chi2inv(0.99999,2)],...
                'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
            
        end
        
        function plotProteinFit(this)
            
            index = find(this.proteinFit.chi == min(this.proteinFit.chi(:)), 1);
            [m, n] = ind2sub(size(this.ed.profiles), index);
            this.plotDataAndFit(this.proteinFit.refnorm, this.proteinFit.ref_fit{m, n});
            
        end
        
    end
    
    methods(Static)
        
        function pmap = pmap()
            
            pmap = [0,0,0.562500000000000;0,0,0.598958313465118;0,0,0.635416686534882;...
                0,0,0.671875000000000;0,0,0.708333313465118;0,0,0.744791686534882;...
                0,0,0.781250000000000;0,0,0.817708313465118;0,0,0.854166686534882;...
                0,0,0.890625000000000;0,0,0.927083313465118;0,0,0.963541686534882;...
                0,0,1;0,0.200000002980232,1;0,0.400000005960465,1;0,0.600000023841858,1;...
                0,0.800000011920929,1;0,1,1;0.142857149243355,1,0.857142865657806;...
                0.285714298486710,1,0.714285731315613;0.428571432828903,1,0.571428596973419;...
                0.571428596973419,1,0.428571432828903;0.714285731315613,1,0.285714298486710;...
                0.857142865657806,1,0.142857149243355;1,1,0;1,0.946708440780640,0;...
                1,0.893416941165924,0;1,0.840125381946564,0;1,0.786833882331848,0;...
                1,0.733542323112488,0;1,0.680250763893127,0;1,0.626959264278412,0;...
                1,0.573667705059052,0;1,0.520376205444336,0;1,0.467084646224976,0;...
                1,0.413793116807938,0;1,0.331034481525421,0;1,0.248275876045227,0;...
                1,0.165517240762711,0;1,0.0827586203813553,0;1,0,0;0.978260874748230,0,0;...
                0.956521749496460,0,0;0.934782624244690,0,0;0.913043498992920,0,0;...
                0.891304373741150,0,0;0.869565188884735,0,0;0.847826063632965,0,0;...
                0.826086938381195,0,0;0.804347813129425,0,0;0.782608687877655,0,0;...
                0.760869562625885,0,0;0.739130437374115,0,0;0.717391312122345,0,0;...
                0.695652186870575,0,0;0.673913061618805,0,0;0.652173936367035,0,0;...
                0.630434811115265,0,0;0.608695626258850,0,0;0.586956501007080,0,0;...
                0.565217375755310,0,0;0.543478250503540,0,0;0.521739125251770,0,0;...
                0.500000000000000,0,0];
            
        end
        
        function plotDataAndFit(ref, refFit)
            
            figure;
            errorbar(ref(:, 1), ref(:, 2), ref(:, 3), '.', 'color', [0, 0, 0] + 0.5, 'linewidth', 1.2);
            hold on;
            plot(refFit(:, 1), refFit(:, 3), '-k', 'linewidth', 2.4);
            set(gca, 'fontsize', 14);
            hold off;
            xlabel('$$ Q_z (\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('Normalized Reflectivity', 'fontsize', 16);
            legend('Data', 'Fit');
            
        end
        
    end
    
end