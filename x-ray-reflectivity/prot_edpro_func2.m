function [ED_prof, R] = prot_edpro_func2(pdbfile, rotlist, pdbflag, chainflag) %, currentFolder, pathname, pdbflag)
%This program is for calculating electron density
%of protein inside a smallest box with air or with buffer
%by using grid method
%Written by Chiuhao Chen and Shekhar Garde  Nov 5,2006
%modified on March/11/2008 by Chiuhao  Chen
%converted to matlab .m file on Aug/10/2012 by Zhiliang Gong
%rewritten for better functionality and to implement in a larger
%programming package Oct/8/2014 by Gregory T. Tietjen
%modified as a function for programming package Jan. 2015 by Daniel Kerr
%Completely revamped from scratch Sept. 2016 by Daniel Kerr

format compact;
fslice  =  0.5;     %fine slice 0.5A in x,y,z directions
    
    % Note that with this textscan call the format is for pdb writen from a
    % selection in VMD: i.e. ATOM # Type Res Seg # x y z # # PROT
    fid = fopen(pdbfile);
    if chainflag == 1
        raw = textscan(fid,'%s %d %s %*s %*s %*d %f %f %f %*f %*f %*s', 'HeaderLines', 1);
    else
        raw = textscan(fid,'%s %d %s %*s %*d %f %f %f %*f %*f %*s', 'HeaderLines', 1);
    end
    fclose(fid);
    %cd(currentFolder)
    
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
        elseif adata.element(k) == 'P'
                adata.radius(k) = 1.95;
                adata.electrons(k) = 15.0;
        elseif adata.element(k) == 'C' && strcmp(tname{k}(1:3), 'CAL') == 1
                adata.radius(k)= 2.4;
                adata.electrons(k) = 20.0;
        else
                sprintf('Unable to identify atom type: atomnumber = %d, name %s ',k ,adata.element(k))
        end
        
        etotpdb = etotpdb + adata.electrons(k);
        
    end

%Now start to set the x,y,z axis
    
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
        theta = theta_rot*pi/180;
        
        % The theta rotation matrix about the x axis following right hand rule is applied next using the matrix:
        % [ 1     0           0      ]
        % [ 0 cos(theta) -sin(theta) ]
        % [ 0 sin(theta)  costheta)  ]
        % 
        
        x = adata.x*cos(phi) - adata.y*sin(phi);
        y = adata.x*cos(theta)*sin(phi) + adata.y*cos(theta)*cos(phi) - adata.z*sin(theta);
        z = adata.x*sin(theta)*sin(phi) + adata.y*sin(theta)*cos(phi) + adata.z*cos(theta);
        
        if pdbflag == 1
            outf = sprintf('t%03dp%03d.pdb',theta_rot,phi_rot);
            wrf2 = fopen(outf,'a+');
            if wrf2  == -1
                sprintf('Cannot create data file: %s',outf)
            end
        
            %cd(pathname)
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
            %cd(currentFolder)
            fclose(wrf2);
        end
        
        
        %Find the box size of protein at certain orientation
        %X; Y; Z  directions
        xtop = max(x+adata.radius);
        xbot = min(x-adata.radius);
        ytop = max(y+adata.radius);
        ybot = min(y-adata.radius);
        ztop = max(z+adata.radius);
        zbot = min(z-adata.radius);
        sprintf('xy area of the box = %f\n ',abs(xtop-xbot)*abs(ytop-ybot))
        
        %Grid array, elements will be # of electrons in each grid. Each
        %grid will be indexed by the coordinate of its center.
        
        xfgridnum = round((xtop-xbot)/fslice);
        yfgridnum = round((ytop-ybot)/fslice);
        zfgridnum = round((ztop-zbot)/fslice);
        
        xfslice = (xtop-xbot)/xfgridnum;
        yfslice = (ytop-ybot)/yfgridnum;
        zfslice = (ztop-zbot)/zfgridnum;
        fvol = xfslice*yfslice*zfslice;
        
        Elec_grid = zeros(xfgridnum,yfgridnum,zfgridnum);
        xf_grid = zeros(xfgridnum,yfgridnum,zfgridnum);
        yf_grid = zeros(xfgridnum,yfgridnum,zfgridnum);
        zf_grid = zeros(xfgridnum,yfgridnum,zfgridnum);
        
        for x_ind = 1:xfgridnum
            for y_ind = 1:yfgridnum
                for z_ind = 1:zfgridnum
                    xf_grid(x_ind,y_ind,z_ind) = xbot + xfslice*(x_ind-0.5);
                    yf_grid(x_ind,y_ind,z_ind) = ybot + yfslice*(y_ind-0.5);
                    zf_grid(x_ind,y_ind,z_ind) = zbot + zfslice*(z_ind-0.5);
                end
            end
        end
        tic
        for atomcount = 1:length(adata.element)
            %tic
            %Array containing the current atom's electrons divided evenly
            %among the grids the atom's van der waal's sphere intersects
            
            Elec_frac = zeros(xfgridnum,yfgridnum,zfgridnum);
            
            %Need to find closest and farthest corners of the grid
            
                  %7%%%%%%%%%%%%%%%%%%8
               %%%%%%%%%%%%%%%%%%%%%  %
%           3%%%%%%%%%%%%%%%%%%%4     %
            %%%%%%%%%%%%%%%%%%%%%     %
            %%%%%%%%%%%%%%%%%%%%%     %
            %%%%%%%%%%%%%%%%%%%%%     %
            %%%%%%%5%%%%%%%%%%%%%%%%%%6
            %%%%%%%%%%%%%%%%%%%%%%%%
%           1%%%%%%%%%%%%%%%%%%%2
            
%             corner(:,:,:,1) = ((x_grid-xslice) - adata.x(atomcount)).^2 + ((y_grid-yslice) - adata.y(atomcount)).^2 +((z_grid-zslice) - adata.z(atomcount)).^2;
%             corner(:,:,:,2) = (x_grid - adata.x(atomcount)).^2 + ((y_grid-yslice) - adata.y(atomcount)).^2 +((z_grid-zslice) - adata.z(atomcount)).^2;
%             corner(:,:,:,3) = ((x_grid-xslice) - adata.x(atomcount)).^2 + (y_grid - adata.y(atomcount)).^2 +((z_grid-zslice) - adata.z(atomcount)).^2;
%             corner(:,:,:,4) = (x_grid - adata.x(atomcount)).^2 + (y_grid - adata.y(atomcount)).^2 +((z_grid-zslice) - adata.z(atomcount)).^2;
%             corner(:,:,:,5) = ((x_grid-xslice) - adata.x(atomcount)).^2 + ((y_grid-yslice) - adata.y(atomcount)).^2 +(z_grid - adata.z(atomcount)).^2;
%             corner(:,:,:,6) = (x_grid - adata.x(atomcount)).^2 + ((y_grid-yslice) - adata.y(atomcount)).^2 +(z_grid - adata.z(atomcount)).^2;
%             corner(:,:,:,7) = ((x_grid-xslice) - adata.x(atomcount)).^2 + (y_grid - adata.y(atomcount)).^2 +(z_grid - adata.z(atomcount)).^2;
%             corner(:,:,:,8) = (x_grid - adata.x(atomcount)).^2 + (y_grid - adata.y(atomcount)).^2 +(z_grid - adata.z(atomcount)).^2;
            
            atomdist = (xf_grid - adata.x(atomcount)).^2 + (yf_grid - adata.y(atomcount)).^2 + (zf_grid - adata.z(atomcount)).^2;
            %Find the grids with at least one corner contained within the
            %atom's van der waal's sphere
            
            ind_list = (atomdist <= adata.radius(atomcount)^2); %&& (max(corner1,corner2,corner3,corner4,corner5,corner6,corner7,corner8) >= adata.radius(atomcount));
            
            Elec_frac(ind_list) = (adata.electrons(atomcount)*fvol)/((4/3)*pi*adata.radius(atomcount)^3);
            
            Elec_grid = Elec_grid + Elec_frac;
            
            clear Elec_frac;
            %toc
        end
        toc
        ED_prof = zeros(zfgridnum,3);
        parray = zeros(zfgridnum,1);
        barray = zeros(zfgridnum,1);
        earray = zeros(zfgridnum,1);
        tote = 0;
        
        numgridarea = any(Elec_grid,3);
        minarea = sum(numgridarea(:));
        
        for lay = 1:zfgridnum
                
            ED_prof(lay,1) = zf_grid(1,1,lay);
            %ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*xfgridnum*yfgridnum);
            ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*minarea);
            ED_prof(lay,3) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol/(fvol*xfgridnum*yfgridnum);

            parray(lay) = nnz(Elec_grid(:,:,lay))*fvol;
            barray(lay) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol;
            earray(lay) = sum(sum(Elec_grid(:,:,lay)));
            tote = tote + earray(lay);
        end
        
        ED_prof = flipud(ED_prof);
        
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
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xfgridnum = xfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).yfgridnum = yfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).zfgridnum = zfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).fvol = fvol;
    end
end

end