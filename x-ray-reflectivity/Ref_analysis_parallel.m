function [ R ] = Ref_analysis_parallel(  )
%Master program for complete reduction and error analysis of x-ray
%reflectivity data.

%Embedded function list:
    % qzoff_fit
    % Rfcalc3
    % Qc_calc
    % Qz_calc
    % lipid2box
    % Lipid_Prot_EDcalc4
    % parratt4
    % Ref_reduce3_man
    % refsave
    % prot_edpro_func_grid
    % chimapplot

%Author Gregory T. Tietjen 01/05/14
%CoAuthor Daniel Kerr 09/22/14
%Version 1.3
%Updated 05/07/14 to include lipid only fit data in output R structure
%Updated 11/17/14 to include error analysis for lipid only fits
%Updated 02/05/15 to include protein analysis
%Updated 04/01/16 completely overhauled removing cplot functionality

%%% Enter While loop for selecting desired functionality %%%

cont = 1; % Continuation variable; this is 1 to continue program and 0 to exit.

opts = optimset('Display', 'off'); %Settings for fitting so fitting iterations are not printed to terminal

%%% Choose initial functionality %%%
option = input('\nWhat would you like to do: \n Enter 1 for qz offset fit \n Enter 2 for Data fitting \n Enter 3 for protein ED profile generation \n=> ');


while cont == 1
     

        %%% Qz offset fitting routine %%%
    if option == 1 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Get model type  %%%
        model_type = input('\nWhat model will we be using? \n Enter 1 for lipid only \n Enter 2 for lipid + protein \n =>');
        %disp('test')

        %%% Step 1: Raw Data import and determine data model type %%%

            %Get raw Reflectivity data
                [FileName,PathName,~]= uigetfile('.txt','Select the raw reflectivity file');
                cd(PathName);
                R.rawdata.ref = importdata(FileName);
                R.input.filename = FileName;
                R.input.filepath = PathName;
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%% Step 2: Determine the qzoffset %%%

            %First determine best 'naive' guess of qc_observed
                %***Note that this will always be an underestimate of the actual 
                %   offset and so a model dependent fitting must be used to 
                %   identify the optimal offset.***

                R.input.p_buff = input('\nInput supbphase electron density value: \n =>');
                R.input.qc = Qc_calc(R.input.p_buff);

                qzoff_naive = qzoff_fit(R.rawdata.ref,R.input.p_buff );

                % Need to round this to 5 decimal places to keep code clean and for
                % file naming purposes

                qzoff_naive = floor(qzoff_naive*10^5)*10^-5;
                
                qc_cutoff_ind = find(R.rawdata.ref(:,1)<=0.026);

            %Use qzoff_naive as a minimum value for determining true
            %qzoffset and use fitting to determine optimal value
            
            %Select fit parameters based on previously chosen model - note that
            %roughness will not be used as a fit parameter at this
            %point to avoid any couple fit parameter issues
            
            if model_type == 1
                x0 = [qzoff_naive, 13, 10.4, 0.223, 0.453]; %qz, l_tail, l_head, p_tail, p_head
                lb = [-Inf, 0, 0, 0, 0];
                ub = [Inf, 100, 30, 100, 100];
                ProtFlag = 0;
                qzoff_edprofile = [];
            elseif model_type == 2
                x0 = [qzoff_naive, -10, 0.35, 13, 10.4, 0.223, 0.453]; %qz, prot_pos, cov, l_tail, l_head, p_tail, p_head
                lb = [-Inf, -100, 0, 0, 0, 0, 0];
                ub = [Inf, 0, 1, 100, 30, 100, 100];
                ProtFlag = 1;
                
                currentFolder = pwd;
                [edfile,PathName,~]  = uigetfile('*.ed','Select the ED file you want to read.');
                if strcmp(PathName, [currentFolder, '\']) == 0
                    cd(PathName);
                    qzoff_edprofile = importdata(edfile);
                    cd(currentFolder);
                else
                    qzoff_edprofile = importdata(edfile);
                end
            end
                
                qzoff_fit_func = @(x, xdata)Qz_calc(x, xdata, qzoff_edprofile, R.input.p_buff, 3.4, R.input.qc, ProtFlag, R.rawdata.ref);
                
                [x, resnorm] = lsqcurvefit(qzoff_fit_func, x0, R.rawdata.ref([2:4 (length(qc_cutoff_ind)+1):end],1),...
                    zeros(size(R.rawdata.ref([2:4 (length(qc_cutoff_ind)+1):end],1))), lb, ub, opts);
                
                qzoff_min = x(1);
                R.offset.chisqrd = resnorm;
                
                [ R.offset.ref, R.offset.refnorm, R.offset.refnormcut,~,~] = ...
                    Ref_reduce3_man(R.rawdata.ref,0.026,qzoff_min,R.input.qc);
                
                R.input.qzoffset = qzoff_min;
                R.offset.ref = R.rawdata.ref;
                R.offset.ref(:,1) = R.offset.ref(:,1) - R.input.qzoffset;
                %[ED, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), R.input.p_buff, 3.4, ProtFlag);
                %R.offset.reflparratt = parratt4(ED, ddlay, R.rawdata.ref([2:4 (length(qc_cutoff_ind)+1):end],1)-x(1), R.input.p_buff);
                fprintf('\n%s %f \n','The Qzoffset is:', R.input.qzoffset)


        filename = input('Give a filename for saving the matfile: \n =>','s');
        save(filename,'R');
        
        
        

        
        %%% Data fitting routine %%%
    elseif option == 2 
            
            
            
        if exist('R')==0
                
            [FileName,PathName,~]= uigetfile('.mat','Select a saved R structure mat file (these must be from the qzoffset feature of this program):');
            cd(PathName);
            R = importdata(FileName);
                
        end
        
        data_name = input('\n What is the name of the dataset? \n=>','s');
            
            %%% Get model type  %%%
            model_type = input('\nWhat model will we be using? \n Enter 1 for lipid only \n Enter 2 for lipid + protein \n =>');

            
            %First prepare data for fitting by removing bad data points 
            %This can be done either by cutting all data below ~0.026 or by
            %selectively removing bad data points (not yet functional)
            
            fprintf('\n %s\n', 'The current version cuts all data below 0.026 for fitting purposes - this data is typically unreliable') 
            cutind = find(R.offset.refnorm(:,1)>0.026,1,'first');
            R.fitdata.refnorm = R.offset.refnorm(cutind+1:end,:);
            
            if model_type == 1
                
                %Lipid only fitting
                    
                    fit_sigma = input('\nDo you want sigma to be a fit parameter? \n (Enter 1 for yes or 0 for no) \n =>');
                    
                    if fit_sigma == 1
                        
                        mkdir box2fit-sigmafree
                        cd box2fit-sigmafree
                        sigma = [];
                        x0 = [3.4, 13, 10.4, 0.223, 0.453]; %sigma, l_tail, l_head, p_tail, p_head
                        lb = [0, 0, 0, 0, 0];
                        ub = [100, 100, 30, 100, 100];
                        
                    elseif fit_sigma == 0 
                        
                        mkdir box2fit-sigmafixed
                        cd box2fit-sigmafixed
                        sigma = input('\nThen what value for sigma? \n =>');
                        x0 = [13, 10.4, 0.223, 0.453]; %l_tail, l_head, p_tail, p_head
                        lb = [0, 0, 0, 0];
                        ub = [100, 30, 100, 100];
                        
                    end
                    
                    lipid2box_func = @(x,xdata)lipid2box(x, xdata, R.fitdata.refnorm(:,3), R.input.p_buff, fit_sigma, sigma);
                    [x, R.fitdata.chi2] = lsqcurvefit(lipid2box_func, x0, R.fitdata.refnorm(:,1), R.fitdata.refnorm(:,2)./R.fitdata.refnorm(:,3), lb, ub, opts);
                    
                    if fit_sigma == 1
                        R.fitdata.l_tail = x(2);
                        R.fitdata.l_head = x(3);
                        R.fitdata.p_tail = x(4);
                        R.fitdata.p_head = x(5);
                        R.fitdata.sigma = x(1);
                    elseif fit_sigma == 0
                        R.fitdata.l_tail = x(1);
                        R.fitdata.l_head = x(2);
                        R.fitdata.p_tail = x(3);
                        R.fitdata.p_head = x(4);
                        R.fitdata.sigma = sigma;
                    end
                    
                    R.fitdata.p_buff = R.input.p_buff;

                    %Finally generate electron density and bet fit
                    %reflectivity                    
                    
                    if fit_sigma == 1
                        [R.fitdata.ed, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), R.fitdata.p_buff, x(1) , 0);
                    elseif fit_sigma == 0
                        [R.fitdata.ed, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(1), x(2), x(3), x(4), R.fitdata.p_buff, R.fitdata.sigma, 0);
                    end
                        
                    R.fitdata.ref_fit = parratt4(R.fitdata.ed, ddlay, R.fitdata.refnorm(:,1), R.fitdata.p_buff);
                    
                    filename = input('Give a filename for saving the matfile: \n =>','s');
                    save(filename,'R');
                    
                    
            elseif model_type == 2
                %Lipid+Protein 2 box fitting
                
                model_name = input('Name of protein structure? \n =>', 's');
                fit_sigma = input('\nDo you want to fit sigma? \n (Enter 1 for yes or 0 for no) \n =>');
                    
                    if fit_sigma == 1
                        
                        x0 = [3.4, -10, 0.35, 13, 10.4, 0.223, 0.453]; %sigma, prot_pos, cov, sigma, l_tail, l_head, p_tail, p_head
                        lb = [0, -100, 0, 0, 0, 0, 0];
                        ub = [100, 0, 1, 100, 30, 100, 100];
                        
                    elseif fit_sigma == 0 
                        
                        sigma = input('\nThen what value for sigma? \n =>');
                        x0 = [-10, 0.35, 13, 10.4, 0.223, 0.453]; %prot_pos, cov, l_tail, l_head, p_tail, p_head
                        lb = [-100, 0, 0, 0, 0, 0];
                        ub = [0, 1, 100, 30, 100, 100];
                    end
                
                %Retrieve ed profile files
                
                originalFolder = pwd;
                [~,PathName,~]  = uigetfile('*.ed','Select a file from the ed profiles folder that will be used in this fit');
                if strcmp(PathName, [originalFolder, '\']) == 0
                    cd(PathName);
                end
                rot = input('Please input the lower bound, step, and higher bound of theta and phi angles as consecutive array: \n =>');
                
                edcheck = zeros(length(rot(1):rot(2):rot(3)), length(rot(4):rot(5):rot(6)));
                theta_rotcount = 0;
                for theta_rot = rot(1):rot(2):rot(3)
                    theta_rotcount = theta_rotcount + 1;
                    phi_rotcount = 0;
                    for phi_rot = rot(4):rot(5):rot(6)
                        phi_rotcount = phi_rotcount + 1;
                        edcheck(theta_rotcount,phi_rotcount) = exist(sprintf('t%03dp%03d.ed', theta_rot, phi_rot), 'file');
%                         if edcheck(theta_rotcount,phi_rotcount) == 1
%                             proted.(sprintf('t%03dp%03d', theta_rot, phi_rot)) = importdata(sprintf('t%03dp%03d.ed', theta_rot, phi_rot));
%                         end
                    end
                end
                
                edchecksize = size(edcheck);
                
                if nnz(edcheck) ~= edchecksize(1)*edchecksize(2)
                    edgen = input('Some ed profiles are missing from the folder \n Enter 1 to generate missing ed profiles \n Enter 2 to skip fitting missing ed profiles \n =>');                    
                elseif nnz(edcheck) == edchecksize(1)*edchecksize(2)
                    edgen = 0;
                end
                
                theta_rotcount = 0;
                thetacount = 0;
                for theta_rot = rot(1):rot(2):rot(3)
                    phi_rotcount = 0;
                    phicount = 0;
                    theta_rotcount = theta_rotcount + 1;
                    if nnz(edcheck(theta_rotcount, :)) ~= length(edcheck(theta_rotcount, :))
                        thetacount = thetacount + 1;
                    end
                    for phi_rot = rot(4):rot(5):rot(6)
                        phi_rotcount = phi_rotcount + 1;
                        if edcheck(theta_rotcount,phi_rotcount) == 0
                            phicount = phicount + 1;
                            rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)]) = [theta_rot, phi_rot];
                        end
                    end
                end
                
                if edgen == 1
                    [pdbfile,PathName,~]  = uigetfile('*.pdb','Select the PDB file you want to read.');
                    cd(PathName);
                    R.edprof = prot_edpro_func_grid(pdbfile, rotlist, currentFolder, PathName, 0);
                    clear rotlist
                    theta_rotcount = 0;
                    for theta_rot = rot(1):rot(2):rot(3)
                        phi_rotcount = 0;
                        theta_rotcount = theta_rotcount + 1;
                        for phi_rot = rot(4):rot(5):rot(6)
                            phi_rotcount = phi_rotcount + 1;
                            rotlist.(['t', num2str(theta_rotcount)]).(['p', num2str(phi_rotcount)]) = [theta_rot, phi_rot];
                        end
                    end
                    
                    cd(currentFolder);
                    
                elseif edgen == 2
                    theta_rotcount = 0;
                    thetacount = 0;
                    for theta_rot = rot(1):rot(2):rot(3)
                        phi_rotcount = 0;
                    phicount = 0;
                        theta_rotcount = theta_rotcount + 1;
                        if nnz(edcheck(theta_rotcount, :)) ~= length(edcheck(theta_rotcount, :))
                            thetacount = thetacount + 1;
                        end
                        for phi_rot = rot(4):rot(5):rot(6)
                            phi_rotcount = phi_rotcount + 1;
                            if edcheck(theta_rotcount,phi_rotcount) == 1
                                phicount = phicount + 1;
                                rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)]) = [theta_rot, phi_rot];
                            end
                        end
                    end
                elseif edgen == 0
                    theta_rotcount = 0;
                    for theta_rot = rot(1):rot(2):rot(3)
                        phi_rotcount = 0;
                        theta_rotcount = theta_rotcount + 1;
                        for phi_rot = rot(4):rot(5):rot(6)
                            phi_rotcount = phi_rotcount + 1;
                            rotlist.(['t', num2str(theta_rotcount)]).(['p', num2str(phi_rotcount)]) = [theta_rot, phi_rot];
                        end
                    end
                end
                    
            %Fitting Procedure
                
            tic;
                refnorm = R.fitdata.refnorm;
                p_buff = R.input.p_buff;
                
                anglelength = 0;
                count = 0;
                
                namestheta = fieldnames(rotlist);
                numtheta = length(namestheta);
                
                for thetacount = 1:numtheta
                    namesphi = fieldnames(rotlist.(['t', num2str(thetacount)]));
                    numphi = length(namesphi);
                    anglelength = anglelength + numphi;
                    for phicount = 1:numphi
                        count = count +1;
                        theta_rot(count) = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(1);
                        phi_rot(count) = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(2);
                    end
                end
                
                x = zeros(anglelength,length(x0));
                chi = zeros(anglelength);
                
                parfor count = 1:anglelength
                    
                    prot_ed = importdata(sprintf('t%03dp%03d.ed', theta_rot(count), phi_rot(count)));
                        
                    prot2box_func = @(x,xdata)prot2box(x, xdata, refnorm(:,3), prot_ed, p_buff, fit_sigma, sigma);
                    [x(count,:), chi(count)] =...
                            lsqcurvefit(prot2box_func, x0, refnorm(:,1), refnorm(:,2)./refnorm(:,3), lb, ub, opts);
                end
                
                count = 0;
                
                for thetacount = 1:numtheta
                    namesphi = fieldnames(rotlist.(['t', num2str(thetacount)]));
                    numphi = length(namesphi);
                    for phicount = 1:numphi
                        theta_rot = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(1);
                        phi_rot = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(2);
                        
                        prot_ed = importdata(sprintf('t%03dp%03d.ed', theta_rot, phi_rot));
                        
                        count = count + 1;

                        if fit_sigma == 1
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).chi = chi(count);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).protpos = x(count,2);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).cov = x(count,3);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_tail = x(count,4);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_head = x(count,5);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_tail = x(count,6);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_head = x(count,7);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_buff = R.input.p_buff;
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).sigma = x(count,1);
                        elseif fit_sigma == 0
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).chi = chi(count);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).protpos = x(count,1);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).cov = x(count,2);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_tail = x(count,3);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_head = x(count,4);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_tail = x(count,5);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_head = x(count,6);
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_buff = R.input.p_buff;
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).sigma = sigma;
                        end

                        %Finally generate electron density and best fit
                        %reflectivity

                        [R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).ed, ddlay] =  Lipid_Prot_EDcalc4(prot_ed, R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).protpos,...
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).cov, R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_tail,...
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).l_head, R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_tail,...
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_head, R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_buff,...
                            R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).sigma, 1);


                        R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).ref_fit = parratt4(R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).ed,...
                            ddlay, R.fitdata.refnorm(:,1), R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).p_buff);

                        clear prot_ed

                    end
                end

                R = chimapplot(R, rot);
                toc;

                cd(originalFolder); 
                
                filename = ['prot2box_', data_name, '_', model_name, '_datafit'];
                customfilename = 0;
                
                while ismember(customfilename, [1,2]) == 0
                    customfilename = input(['Default file name is "', filename, '" \n Enter 1 to accept filename \n Enter 2 to enter custom name \n =>']);
                end
                
                if customfilename == 2
                    filename = input('Give a filename for saving the matfile: \n =>','s');
                end
                save(filename,'R');
                
            end
                                
    elseif option == 3
        [pdbfile,PathName,~]  = uigetfile('*.pdb','Select the PDB file you want to read.');
        currentFolder = pwd;
        if isempty(strfind(currentFolder, strtok(pdbfile, '.'))) == 1
            mkdir(strtok(pdbfile, '.'));
            cd(strtok(pdbfile, '.'));
            currentFolder = pwd;
        end
        rot = input('Please input the lower bound, step, and higher bound of theta and phi angles as consecutive array: \n =>');
        pdbflag = input('Produce pdb files for each orientation? \n 1 Yes \n 2 No \n =>');
        edcheck = zeros(length(rot(1):rot(2):rot(3)), length(rot(4):rot(5):rot(6)));
        theta_rotcount = 0;
        for theta_rot = rot(1):rot(2):rot(3)
            theta_rotcount = theta_rotcount + 1;
            phi_rotcount = 0;
            for phi_rot = rot(4):rot(5):rot(6)
                phi_rotcount = phi_rotcount + 1;
                edcheck(theta_rotcount,phi_rotcount) = exist(sprintf('t%03dp%03d.ed', theta_rot, phi_rot), 'file');
            end
        end
        
        if nnz(edcheck) > 0
            overwrite = input('Overwrite existing ed profiles? \n Enter 1 to overwrite \n Enter 0 to skip files \n =>');
        else
            overwrite = 1;
        end
        
        theta_rotcount = 0;
        thetacount = 0;
        for theta_rot = rot(1):rot(2):rot(3)
            phi_rotcount = 0;
            phicount = 0;
            theta_rotcount = theta_rotcount + 1;
            if nnz(edcheck(theta_rotcount, :)) ~= length(edcheck(theta_rotcount, :)) || overwrite == 1
                thetacount = thetacount + 1;
            end
            for phi_rot = rot(4):rot(5):rot(6)
                phi_rotcount = phi_rotcount + 1;
                if edcheck(theta_rotcount,phi_rotcount) == 0 || overwrite == 1
                    phicount = phicount + 1;
                    rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)]) = [theta_rot, phi_rot];
                    if overwrite == 1 && edcheck(theta_rotcount,phi_rotcount) == 1
                        delete(sprintf('t%03dp0%3d.ed', theta_rot, phi_rot));
                    end
                end
            end
        end
        edchecksize = size(edcheck);
        
        if nnz(edcheck) ~= edchecksize(1)*edchecksize(2) || overwrite == 1
            %if overwrite == 1
            %    delete('t*.ed');
            %end
            cd(PathName);
            R.edprof = prot_edpro_func_grid(pdbfile, rotlist, currentFolder, PathName, pdbflag);
        elseif nnz(edcheck) == edchecksize(1)*edchecksize(2)
            cd(currentFolder);
        end
                
    end
    

%%% Either exit program or repeat        
option = input('What would you like to do: \n Enter 1 for qz offset fit \n Enter 2 for Data fitting \n Enter 3 for protein ED profile generation \n Enter 4 to exit the program \n => ');        
      
    if option ~= 4
        
        cont = 1; 
        
    else
        
        cont = 0;
        
    end
    
end   
    

        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ qzoff ] = qzoff_fit(Ref,p_buff )
%Function to fit qz offset of Reflectivity data
%Author: Gregory T Tietjen 10/08/2012

%Calculate qc

qc = Qc_calc(p_buff);

%Determine data points below qc:

ind = find(Ref(:,2)>.95*Ref(1,2));
norm = mean(Ref(ind,2));

qz = Ref(:,1);
refdata = Ref(:,2)/norm;

%Define fit function
%ft = fittype('Rfcalc3(x,qc,qzoff)','problem','qc');
Rfcalcfit = @(qzoff,xdata)Rfcalc3(xdata,qc,qzoff);

%Perform fit and extract qzoff
%f = fit(qz,refdata,ft, 'problem',qc,'StartPoint', .0004);
%qzoff = coeffvalues(f);
opts = optimset('Display', 'off');
qzoff = lsqcurvefit(Rfcalcfit, 0.0004, qz, refdata, [], [], opts);
    
%plot(f,qz,refdata)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [ oRef, oRefoRf, ocRefoRf, qc, qzoff] = Ref_reduce3_man(Ref, qc_cut, qzoff, qc)
%Function for processing Raw Reflectivity Data



format long

%Determine data points below qc:

ind = find(Ref(:,1)<=qc_cut);

    
%Rf generation

   
    Rf = Rfcalc3(Ref(:,1),qc,qzoff);
    
    oRef = Ref;
    oRef(:,1) = oRef(:,1)-qzoff;
    oRefoRf(:,1) = Ref(:,1)-qzoff;
    oRefoRf(:,2) = Ref(:,2)./Rf;
    oRefoRf(:,3) = Ref(:,3)./Rf;
        
  
%Rf cut generation
    
    ocRefoRf = oRefoRf([2:4 (length(ind)+1):end],:);      

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R] = Qz_calc(x, reflxaxis, prot_ed, p_buff, sigma, qc, ProtFlag, refl)

%Function for fitting qz offset to 2 box model with/without protein

[~, ~, Refl_Fresnel_qz_shift_cutoff, ~, ~] = Ref_reduce3_man(refl, 0.026, x(1), qc);

if ProtFlag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, sigma, ProtFlag);
elseif ProtFlag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4([], [], [], x(2), x(3), x(4), x(5), p_buff, sigma, ProtFlag);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis-x(1), p_buff);

R = (Refl_Fresnel_qz_shift_cutoff(:,2) - Reflparratt(:,3))./Refl_Fresnel_qz_shift_cutoff(:,3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R] = prot2box(x, reflxaxis, reflerr, prot_ed, p_buff, sigmaflag, sigma)

%Function for fitting lipid 2box model to qz offset reflectivity data

if sigmaflag == 1
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(2), x(3), x(4), x(5), x(6), x(7), p_buff, x(1) , 1);
elseif sigmaflag == 0
    [ED, ddlay] = Lipid_Prot_EDcalc4(prot_ed, x(1), x(2), x(3), x(4), x(5), x(6), p_buff, sigma, 1);
end

Reflparratt = parratt4(ED, ddlay, reflxaxis, p_buff);

R = Reflparratt(:,3)./reflerr;

end

function [R] = prot_edpro_func_grid(pdbfile, rotlist, currentFolder, pathname, pdbflag)
%This program is for calculating electron density
%of protein inside a smallest box with air or with buffer
%by using grid method
%Written by Chiuhao Chen and Shekhar Garde  Nov 5,2006
%modified on March/11/2008 by Chiuhao  Chen
%converted to matlab .m file on Aug/10/2012 by Zhiliang Gong
%rewritten for better functionality and to implement in a larger
%programming package Oct/8/2014 by Gregory T. Tietjen
%modified as a function for programming package Jan. 2015 by Daniel Kerr

format compact;
maxatm  =  10000;   %assume total number of atoms
cslice  =  2.0;     %fixed coarse slice 2A in x,y,z directions
fslice  =  0.5;     %fine slice 0.5A in x,y,z directions
cornum  =  100;     %coarse pieces
    
    % Note that with this textscan call the format is for pdb writen from a
    % selection in VMD: i.e. ATOM # Type Res Seg # x y z # # PROT
    fid = fopen(pdbfile);
    raw = textscan(fid,'%s %d %s %*s %*s %*d %f %f %f %*f %*f %*s', 'HeaderLines', 1);
    fclose(fid);
    cd(currentFolder)
    
    %Extract data from raw structure
    vn = raw{1}; %Grab Atom or Anisou label
    sn = raw{2}; %This grabs the atom number, used to distinguish atom from anisou.
    tn = raw{3}; %atom type in list form - note this is in MD parameter file notation; 
    %need to get first letter of each entry
    adata.xraw = raw{4};
    adata.yraw = raw{5};
    adata.zraw = raw{6};
    count = 0;
    for j = 1:length(raw{3})
        if j > 1
            if (strcmp(vn(j), 'ATOM') == 1 ||strcmp(vn(j), 'HETATM') == 1) && sn(j) ~= 0
            %if sn(j) ~= sn(j-1) && sn(j) ~= 0 %Atom and Anisou are labeled by the same number, by assuring that sn(j) ~= sn(j-1), we only grab the atom data
                count = count+1;
                tname{count} = char(tn(j));
                adata.element(count) = tname{count}(1);
                adata.x(count) = adata.xraw(j);
                adata.y(count) = adata.yraw(j);
                adata.z(count) = adata.zraw(j);
            end
        elseif j == 1
            count = 1;
            tname{count} = char(tn(j));
            adata.element(count) = tname{count}(1);
            adata.x(count) = adata.xraw(j);
            adata.y(count) = adata.yraw(j);
            adata.z(count) = adata.zraw(j);
        end
    end
    

    %Next determine atom type and assign radius and # of electrons 
    etotpdb = 0;

    for k  =  1:length(adata.element)
       
        if adata.element(k) == 'C'
                adata.radius(k) = 1.7;
                adata.electrons(k) = 6.0;        
        elseif adata.element(k) == 'N' || adata.element(k) == 'A'
                adata.radius(k)  =  1.55;
                adata.electrons(k)  =  7.0;
        elseif adata.element(k) == 'O'
                adata.radius(k)  =  1.52;
                adata.electrons(k)  =  8.0;
        elseif adata.element(k) == 'S'
                adata.radius(k)  =  1.8;
                adata.electrons(k)  =  16.0;
        elseif adata.element(k) == 'H'
                adata.radius(k)  =  1.1;
                adata.electrons(k)  =  1.0;
        elseif adata.element(k) == 'C' && strcmp(tname{k}(1:3), 'CAL') == 1;
                adata.radius(k)= 2.4;
                adata.electrons(k) = 20.0;
        else
                sprintf('Unable to identify atom type: atomnumber = %d, name %s ',k ,adata.element(k))
        end
        
        etotpdb = etotpdb + adata.electrons(k);
        
    end

%Now start to set the x,y,z axis
x = adata.x;
y = adata.y;
z = adata.z;
    
%Rotate protein with theta and phi angles both between 0 to 2pi, where phi
%is the rotation angle about z axis and theta is the rotation angle about x
%axis

%rot = input('Please input the lower bound, step, and higher bound of theta and phi angles as consecutive array: \n =>');

namestheta = fieldnames(rotlist);
numtheta = length(namestheta);
for thetacount = 1:numtheta
    namesphi = fieldnames(rotlist.(['t', num2str(thetacount)]));
    numphi = length(namesphi);
    for phicount = 1:numphi
        theta_rot = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(1);
        phi_rot = rotlist.(['t', num2str(thetacount)]).(['p', num2str(phicount)])(2);
        
        %theta_rot is the theta angle relative to the normal (rotated about the x axis)
        %phi_rot is the phi rotation about the initial z axis \
        %Note that this version expect the input protein structure to be aligned along its
        %principle axis
        
        % The phi rotation matrix about the Z axis following right hand rule is applied first using the martix:
        % [ cos(phi) -sin(phi) 0 ]
        % [ sin(phi)  cos(phi) 0 ]
        % [ 0         0        1 ]
       
        phi = phi_rot*pi/180;
        x = adata.x*cos(phi) - adata.y*sin(phi);
        y = adata.x*sin(phi) + adata.y*cos(phi);
        z = adata.z;
        %x = xfix;
        %y = yfix;
        %z = zfix;
        
        % The theta rotation matrix about the x axis following right hand rule is applied next using the matrix:
        % [ 1     0           0      ]
        % [ 0 cos(theta) -sin(theta) ]
        % [ 0 sin(theta)  costheta)  ]
        % 
        
        theta = theta_rot*pi/180;
        x = x;
        y = y*cos(theta) - z*sin(theta);
        z = y*sin(theta) + z*cos(theta);
        
        
        xfake = x;
        y = adata.x*cos(theta)*sin(phi) + adata.y*cos(theta)*cos(phi) - adata.z*sin(theta);
        z = adata.x*sin(theta)*sin(phi) + adata.y*sin(theta)*cos(phi) + adata.z*cos(theta);
        
        if pdbflag == 1
            outf = sprintf('t%03dp%03d.pdb',theta_rot,phi_rot);
            wrf2 = fopen(outf,'a+');
            if wrf2  == -1
                sprintf('Cannot create data file: %s',outf)
            end
        
            cd(pathname)
            red = fopen(pdbfile,'r');
            if red  == -1
                sprintf('ERROR: cannot create data file: %s',outf)
                break;
            end
        
        %create the new pdbfile corresponding to the rot1 and rot2 values
            k = 0;
            line = fgets(red, 100);
            while ~strcmp(line,' ') && length(line) > 1
                fir = line(1:4);
                i = strcmp(fir,'ATOM');
                j = strcmp(fir,'HETA');
                if i == 1 || j == 1
                    k = k + 1;
                    befxyz = line(12:30); %the line content before x, y and z coordinates
                    aftxyz = line(55:length(line)); %the line content after x, y and z coordinates
                    fprintf(wrf2,'%s%5d%s%8.3f%8.3f%8.3f%s',[fir, '  '],k,befxyz,x(k),y(k),z(k),aftxyz);
                end
                line = fgets(red, 100);
            end
            fclose(red);
            cd(currentFolder)
            fclose(wrf2);
        end
        
        outf = sprintf('t%03dp%03d.ed',theta_rot,phi_rot);
        wrf1 = fopen(outf,'a+');
        if wrf1  == -1
            sprintf('Cannot create data file: %s',outf)
        end
        
        %Find the box size of protein at certain orientation
        %X; Y; Z  directions
        xtop = max(x+adata.radius);
        maxxatom = find((x+adata.radius)==xtop);
        xbot = min(x-adata.radius);
        minxatom = find((x-adata.radius)==xbot);
        ytop = max(y+adata.radius);
        maxyatom = find((y+adata.radius)==ytop);
        ybot = min(y-adata.radius);
        minyatom = find((y-adata.radius)==ybot);
        ztop = max(z+adata.radius);
        maxzatom = find((z+adata.radius)==ztop);
        zbot = min(z-adata.radius);
        minzatom = find((z-adata.radius)==zbot);
        sprintf('xy area of the box = %f\n ',abs(xtop-xbot)*abs(ytop-ybot))
        
        %Coarse Grid x,y,z directions
        xcgrid = round((xtop-xbot)/cslice);
        ycgrid = round((ytop-ybot)/cslice);
        zcgrid = round((ztop-zbot)/cslice);
        dxc = (xtop-xbot)/xcgrid;
        dyc = (ytop-ybot)/ycgrid;
        dzc = (ztop-zbot)/zcgrid;
        cbox = zeros(xcgrid+1, ycgrid+1, zcgrid+1, cornum);
        elementcheck = zeros(size(adata.radius));
        for k = 1:zcgrid
            for i = 1:xcgrid
                for j = 1:ycgrid
                    zctop = ztop-(k-1)*dzc;
                    zcbot = ztop-k*dzc;
                    yctop = ybot+j*dyc;
                    ycbot = ybot+(j-1)*dyc;
                    xctop = xbot+i*dxc;
                    xcbot = xbot+(i-1)*dxc;
                    t = 0;
                    for m = 1:length(adata.radius)
                        if x(m) >= xcbot && x(m) <= xctop && y(m) >= ycbot && y(m) <= yctop && z(m) >= zcbot && z(m)<= zctop && elementcheck(m) == 0
                            t = t+1;
                            cbox(i,j,k,t) = m;
                            elementcheck(m) = 1;
                        end
                    end
            %if t > cornum
            %    sprintf('too many atoms %d',t)
            %end
                end
            end
        end
        
        %Fine Grid x,y,z directions
        xfgrid = round(dxc/fslice);
        yfgrid = round(dyc/fslice);
        zfgrid = round(dzc/fslice);
        dxf = dxc/xfgrid;
        dyf = dyc/yfgrid;
        dzf = dzc/zfgrid;
        fvol = dxf*dyf*dzf;
        tote = 0.0;
        
        parray = zeros(zcgrid,zfgrid);
        barray = zeros(zcgrid,zfgrid);
        earray = zeros(zcgrid,zfgrid);

        for k = 1:zcgrid
            provol = zeros(1,zfgrid);
            bufvol = provol;
            etot   = provol;
            zhigh  = provol;

            lay    = 0;
            for i = 1:xcgrid
                for j = 1:ycgrid
                    zctop = ztop-(k-1.0)*dzc;
                    zcbot = ztop-k*dzc;
                    yctop = ybot+j*dyc;
                    ycbot = ybot+(j-1.0)*dyc;
                    xctop = xbot+i*dxc;
                    xcbot = xbot+(i-1.0)*dxc;
                    for lay = 1:zfgrid
                        for v = 1:yfgrid
                            for r = 1:xfgrid
                                if j==1 && i==1 && v == 1 && r == 1
                                    zhigh(lay) = zctop - lay*dzf;
                                end
                                zfbot = zctop - lay*dzf;
                                yftop = ycbot + v*dyf;
                                xftop = xcbot + r*dxf;
                                et = 0.0;                           
                                for k1 = (k-1):(k+1)
                                    k1old = k1;
                                    if k1 == k-1 && k == 1
                                        k1 = k;
                                    end
                                    for i1 = (i-1):(i+1)
                                        i1old = i1;
                                        if i1 == i-1 && i == 1
                                            i1 = i;
                                        end
                                        for j1 = (j-1):(j+1)
                                            j1old = j1;
                                            if j1 == j-1 && j == 1
                                                j1 = j;
                                            end
                                            if nnz(cbox(i1,j1,k1,:)) > 0
                                                for h = 1:nnz(cbox(i1,j1,k1,:))
                                                    m = cbox(i1,j1,k1,h);
                                                    if (z(m)-zfbot)^2 + (y(m)-yftop)^2 + (x(m)-xftop)^2 < adata.radius(m)^2
                                                        et = et + (adata.electrons(m)*fvol)/((4/3)*pi*adata.radius(m)^3);
                                                    end
                                                end
                                            end
                                            j1 = j1old;
                                        end
                                        i1 = i1old;
                                    end
                                    k1 = k1old;
                                end
                                if et > 0
                                    provol(lay) = provol(lay)+fvol;
                                    etot(lay) = etot(lay)+et;
                                else
                                    bufvol(lay) = bufvol(lay)+fvol;
                                end
                            end
                        end
                    end
                end
            end
            for lay = 1:zfgrid
                fprintf(wrf1,'%.8f %.8f %.8f\n',zhigh(lay),etot(lay)/(provol(lay)+bufvol(lay)),bufvol(lay)/(provol(lay)+bufvol(lay)));
                tote = tote+etot(lay);
                
                parray(k,lay) = provol(lay);
                barray(k,lay) = bufvol(lay);
                earray(k,lay) = etot(lay);
            end
        end
        fclose('all');
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xyarea = abs(xtop-xbot)*abs(ytop-ybot);
        R.etotpdb = etotpdb;
        display(etotpdb);
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).etot = earray;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).tote = tote;
        display(sum(tote));
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).elecdif = (tote - etotpdb);
        display(sum(tote) - etotpdb);
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).provol = parray;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).bufvol = barray;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ R ] = chimapplot(R, rot)
%Function to create plot of chi^2 values of R structure output of
%Ref_analysis of protein orientations

theta_rotcount = 0;
theta = zeros((rot(3)-rot(1))/rot(2) + 1, 1);
phi = zeros((rot(6)-rot(4))/rot(5) + 1, 1);
chi = zeros(((rot(6)-rot(4))/rot(5)) + 1, (rot(3)-rot(1))/rot(2) + 1);

for theta_rot = rot(1):rot(2):rot(3)
    theta_rotcount = theta_rotcount + 1;
    phi_rotcount = 0;
    theta(theta_rotcount) = theta_rot;
    for phi_rot = rot(4):rot(5):rot(6)
        phi_rotcount = phi_rotcount + 1;
        chi(phi_rotcount, theta_rotcount) = R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).chi;
        if theta_rot == rot(1)
            phi(phi_rotcount) = phi_rot; %only compute phi array once
        end
    end
end

[Chimin, ~] = min(chi(:));

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

dof = length(R.fitdata.refnorm(:,1))-6-1;

R.chimap.theta = theta;
R.chimap.phi = phi;
R.chimap.chiarray = chi;
R.chimap.chinormarray = chi/Chimin;
chisym = sym(R.chimap.chiarray);
R.chimap.likelihoodarray = vpa(exp(-(1/2)*chisym));
R.chimap.likelihoodarray = R.chimap.likelihoodarray/sum(R.chimap.likelihoodarray(:));
R.chimap.likelihoodarray = double(R.chimap.likelihoodarray);
chipvalue = fcdf(R.chimap.chinormarray, dof, dof);
R.chimap.pvaluearray = ones(size(chipvalue)) - chipvalue;
[Iphimin, Ithetamin] = find(chi == Chimin);
if length(Ithetamin) == 1
    R.chimap.min = R.fitdata.(sprintf('t%03dp%03d', theta(Ithetamin), phi(Iphimin)));
    R.chimap.min.theta = theta(Ithetamin);
    R.chimap.min.phi = phi(Iphimin);
else
    for i = 1:length(Ithetamin)
        R.chimap.(['min', num2str(i)]) = R.fitdata.(sprintf('t%03dp%03d', theta(Ithetamin(i)), phi(Iphimin(i))));
        R.chimap.(['min', num2str(i)]).theta = theta(Ithetamin(i));
        R.chimap.(['min', num2str(i)]).phi = phi(Iphimin(i));
    end
end

CorrFact = length(R.chimap.theta)*length(R.chimap.phi);

figure(1)
imagesc(R.chimap.theta,R.chimap.phi,R.chimap.chinormarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap); 
caxis([1 3]); 
colorbar('YTick',[finv(0.95,dof,dof), finv(0.99,dof,dof), finv(0.999,dof,dof), finv(0.9999,dof,dof), finv(0.99999,dof,dof)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by F test P Value', 'FontSize', 20)

figure(2)
imagesc(R.chimap.theta,R.chimap.phi,R.chimap.chinormarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap); 
caxis([1 7])
colorbar('YTick',[finv((1-0.05)^(1/(CorrFact)),dof,dof), finv((1-0.01)^(1/(CorrFact)),dof,dof), finv((1-0.001)^(1/(CorrFact)),dof,dof), finv((1-0.0001)^(1/(CorrFact)),dof,dof), finv((1-0.00001)^(1/(CorrFact)),dof,dof)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by F test P Value Corrected', 'FontSize', 20)

figure(3)
imagesc(R.chimap.theta,R.chimap.phi,(R.chimap.chiarray - Chimin*ones(size(R.chimap.chiarray))));
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap);
caxis([0, chi2inv(0.999999,2)]);
colorbar('YTick',[chi2inv(0.95,2), chi2inv(0.99,2), chi2inv(0.999,2), chi2inv(0.9999,2), chi2inv(0.99999,2)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by \chi^2 P Value', 'FontSize', 20)

figure(4)
surf(R.chimap.theta,R.chimap.phi,R.chimap.likelihoodarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(flipud(pmap));
title('Likelihood Surface Plot', 'FontSize', 20);

end

