function [R] = prot_edpro_func2_covmat(pdbfile, covmat, rotlist, pdbflag, chainflag) %, currentFolder, pathname, pdbflag)
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
        
        rotmat = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)]*[cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1];
        covmat_rot = zeros(3,3,length(adata.element));
        for atomind = 1:length(adata.element)
            covmat_rot(:,:,atomind) = rotmat*covmat(:,:,atomind)*(rotmat^-1);
        end
        
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
        
        outf = sprintf('t%03dp%03d.ed',theta_rot,phi_rot);
        wrf1 = fopen(outf,'a+');
        
        
        %Find the box size of protein at certain orientation
        %X; Y; Z  directions
        xtop = max(x+sqrt(squeeze(covmat_rot(1,1,:)))'+adata.radius) + fslice;
        xbot = min(x-sqrt(squeeze(covmat_rot(1,1,:)))'-adata.radius) - fslice;
        ytop = max(y+sqrt(squeeze(covmat_rot(2,2,:)))'+adata.radius) + fslice;
        ybot = min(y-sqrt(squeeze(covmat_rot(2,2,:)))'-adata.radius) - fslice;
        ztop = max(z+sqrt(squeeze(covmat_rot(3,3,:)))'+adata.radius) + fslice;
        zbot = min(z-sqrt(squeeze(covmat_rot(3,3,:)))'-adata.radius) - fslice;
        sprintf('xy area of the box = %f',abs(xtop-xbot)*abs(ytop-ybot))
        
        xtop_nocov = max(x+adata.radius);
        xbot_nocov = min(x-adata.radius);
        ytop_nocov = max(y+adata.radius);
        ybot_nocov = min(y-adata.radius);
        ztop_nocov = max(z+adata.radius);
        zbot_nocov = min(z-adata.radius);
        area_nocov = abs(xtop_nocov-xbot_nocov)*abs(ytop_nocov-ybot_nocov);
        
        %Grid array, elements will be # of electrons in each grid. Each
        %grid will be indexed by the coordinate of its center.
        
        xfgridnum = round((xtop-xbot)/fslice);
        yfgridnum = round((ytop-ybot)/fslice);
        zfgridnum = round((ztop-zbot)/fslice);
        
        xfgridnum_nocov = round((xtop_nocov-xbot_nocov)/fslice);
        yfgridnum_nocov = round((ytop_nocov-ybot_nocov)/fslice);
        zfgridnum_nocov = round((ztop_nocov-zbot_nocov)/fslice);
        
        xfslice = (xtop-xbot)/xfgridnum;
        yfslice = (ytop-ybot)/yfgridnum;
        zfslice = (ztop-zbot)/zfgridnum;
        fvol = xfslice*yfslice*zfslice;
        
        xfslice_nocov = (xtop_nocov-xbot_nocov)/xfgridnum_nocov;
        yfslice_nocov = (ytop_nocov-ybot_nocov)/yfgridnum_nocov;
        zfslice_nocov = (ztop_nocov-zbot_nocov)/zfgridnum_nocov;
        fvol_nocov = xfslice_nocov*yfslice_nocov*zfslice_nocov;
        
        x_grid_pos = xbot + xfslice * ((1:xfgridnum) - 0.5);
        y_grid_pos = ybot + yfslice * ((1:yfgridnum) - 0.5);
        z_grid_pos = zbot + zfslice * ((1:zfgridnum) - 0.5);
        
        x_diff = repmat(x_grid_pos, length(x), 1) - repmat(x', 1, length(x_grid_pos));
        y_diff = repmat(y_grid_pos, length(y), 1) - repmat(y', 1, length(y_grid_pos));
        z_diff = repmat(z_grid_pos, length(z), 1) - repmat(z', 1, length(z_grid_pos));
        
        [x_top_atom_indices(:,1), x_top_atom_indices(:,2)] = ind2sub([length(x), length(x_grid_pos)], find((x_diff > adata.radius' + sqrt(squeeze(covmat_rot(1,1,:)))) & (x_diff < adata.radius' + sqrt(squeeze(covmat_rot(1,1,:))) + xfslice)));
        [x_bot_atom_indices(:,1), x_bot_atom_indices(:,2)] = ind2sub([length(x), length(x_grid_pos)], find((x_diff < -adata.radius' - sqrt(squeeze(covmat_rot(1,1,:)))) & (x_diff > -adata.radius' - sqrt(squeeze(covmat_rot(1,1,:))) - xfslice)));
        [y_top_atom_indices(:,1), y_top_atom_indices(:,2)] = ind2sub([length(y), length(y_grid_pos)], find((y_diff > adata.radius' + sqrt(squeeze(covmat_rot(2,2,:)))) & (y_diff < adata.radius' + sqrt(squeeze(covmat_rot(2,2,:))) + yfslice)));
        [y_bot_atom_indices(:,1), y_bot_atom_indices(:,2)] = ind2sub([length(y), length(y_grid_pos)], find((y_diff < -adata.radius' - sqrt(squeeze(covmat_rot(2,2,:)))) & (y_diff > -adata.radius' - sqrt(squeeze(covmat_rot(2,2,:))) - yfslice)));
        [z_top_atom_indices(:,1), z_top_atom_indices(:,2)] = ind2sub([length(z), length(z_grid_pos)], find((z_diff > adata.radius' + sqrt(squeeze(covmat_rot(3,3,:)))) & (z_diff < adata.radius' + sqrt(squeeze(covmat_rot(3,3,:))) + zfslice)));
        [z_bot_atom_indices(:,1), z_bot_atom_indices(:,2)] = ind2sub([length(z), length(z_grid_pos)], find((z_diff < -adata.radius' - sqrt(squeeze(covmat_rot(3,3,:)))) & (z_diff > -adata.radius' - sqrt(squeeze(covmat_rot(3,3,:))) - zfslice)));

        x_top_atom_indices = sortrows(x_top_atom_indices);
        x_bot_atom_indices = sortrows(x_bot_atom_indices);
        y_top_atom_indices = sortrows(y_top_atom_indices);
        y_bot_atom_indices = sortrows(y_bot_atom_indices);
        z_top_atom_indices = sortrows(z_top_atom_indices);
        z_bot_atom_indices = sortrows(z_bot_atom_indices);

        x_top_indices = x_top_atom_indices(:,2);
        x_bot_indices = x_bot_atom_indices(:,2);
        y_top_indices = y_top_atom_indices(:,2);
        y_bot_indices = y_bot_atom_indices(:,2);
        z_top_indices = z_top_atom_indices(:,2);
        z_bot_indices = z_bot_atom_indices(:,2);
        
        [xf_grid, yf_grid, zf_grid] = meshgrid(x_grid_pos, y_grid_pos, z_grid_pos);
                    
        xf_grid = permute(xf_grid, [2,1,3]);
        yf_grid = permute(yf_grid, [2,1,3]);
        zf_grid = permute(zf_grid, [2,1,3]);
        
        Elec_grid = zeros(xfgridnum,yfgridnum,zfgridnum);        
        Buffer_grid = ones(xfgridnum,yfgridnum,zfgridnum);
%         atomdistvec = zeros(3,xfgridnum,yfgridnum,zfgridnum);
%         Ellipatomdistvec = zeros(3,xfgridnum,yfgridnum,zfgridnum);
        
        tic
        for atomcount = 1:length(adata.element)
            %tic
            %Array containing the current atom's electrons divided evenly
            %among the grids the atom's van der waal's sphere intersects
            
            %Elec_frac = zeros(xfgridnum,yfgridnum,zfgridnum);
            Elec_frac = zeros(x_top_indices(atomcount)-x_bot_indices(atomcount)+1, y_top_indices(atomcount)-y_bot_indices(atomcount)+1, z_top_indices(atomcount)-z_bot_indices(atomcount)+1);
            Buffer_frac = ones(x_top_indices(atomcount)-x_bot_indices(atomcount)+1, y_top_indices(atomcount)-y_bot_indices(atomcount)+1, z_top_indices(atomcount)-z_bot_indices(atomcount)+1);
            
            if det(covmat_rot(:,:,atomcount)) == 0
                
                atomdist = (xf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - x(atomcount)).^2 ...
                    + (yf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - y(atomcount)).^2 ...
                    + (zf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - z(atomcount)).^2;
                %atomdist = (xf_grid - x(atomcount)).^2 + (yf_grid - y(atomcount)).^2 + (zf_grid - z(atomcount)).^2;
                %Find the grids with centers contained within the
                %atom's van der waal's sphere
                %ind_list = (atomdist <= adata.radius(atomcount)^2);
                %Elec_frac(ind_list) = (adata.electrons(atomcount)*fvol)/((4/3)*pi*adata.radius(atomcount)^3);

                ind_list = (atomdist <= adata.radius(atomcount)^2);

                Elec_frac(ind_list) = (adata.electron(atomcount)*fvol)/((4/3)*pi*adata.radius(atomcount)^3);
                
            else
                [eigvec, eigval] = eig(covmat_rot(:,:,atomcount));
                radiusvec = sqrt(diag(eigval)) + adata.radius(atomcount)*ones(3,1);
                %Add radius of the atom in since it's the atom centers
                %which the covariance matrices represent.
                
                atomdistvec = zeros(3,x_top_indices(atomcount)-x_bot_indices(atomcount)+1, y_top_indices(atomcount)-y_bot_indices(atomcount)+1, z_top_indices(atomcount)-z_bot_indices(atomcount)+1);
                Ellipatomdistvec = zeros(3,x_top_indices(atomcount)-x_bot_indices(atomcount)+1, y_top_indices(atomcount)-y_bot_indices(atomcount)+1, z_top_indices(atomcount)-z_bot_indices(atomcount)+1);
                
                Ellip = eigvec*diag(1./(radiusvec).^2)*(eigvec');
                atomdistvec(1,:,:,:) = (xf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - x(atomcount));
                atomdistvec(2,:,:,:) = (yf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - y(atomcount)); 
                atomdistvec(3,:,:,:) = (zf_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) - z(atomcount));
                %atomdist = atomdistvec.*(Ellip*atomdistvec);
                Ellipatomdistvec(1,:,:,:) = Ellip(1,1)*atomdistvec(1,:,:,:) + Ellip(1,2)*atomdistvec(2,:,:,:) + Ellip(1,3)*atomdistvec(3,:,:,:);
                Ellipatomdistvec(2,:,:,:) = Ellip(2,1)*atomdistvec(1,:,:,:) + Ellip(2,2)*atomdistvec(2,:,:,:) + Ellip(2,3)*atomdistvec(3,:,:,:);
                Ellipatomdistvec(3,:,:,:) = Ellip(3,1)*atomdistvec(1,:,:,:) + Ellip(3,2)*atomdistvec(2,:,:,:) + Ellip(3,3)*atomdistvec(3,:,:,:);
                atomdist = atomdistvec(1,:,:,:).*Ellipatomdistvec(1,:,:,:) + atomdistvec(2,:,:,:).*Ellipatomdistvec(2,:,:,:) +  atomdistvec(3,:,:,:).*Ellipatomdistvec(3,:,:,:);
                %Find the grids with centers contained within the atom's
                %van der waal's sphere smeared over an ellipsoid
                %representing the protein's covariance over the MD
                %trajectory.
                
                ind_list = squeeze((atomdist <= 1));
                
                EllipVol = (4/3)*pi*(sqrt(eigval(1,1))+adata.radius(atomcount))*(sqrt(eigval(2,2))+adata.radius(atomcount))*(sqrt(eigval(3,3))+adata.radius(atomcount));
                VDWVol = (4/3)*pi*adata.radius(atomcount)^3;
                VolRatio = VDWVol/EllipVol;
                
                Elec_frac(ind_list) = (adata.electrons(atomcount)*fvol)/EllipVol;
                %Divide by the volume of the van der waal's smeared
                %covariance ellipse instead of the atom's van der waal's
                %volume since the electrons are smeared over the
                %ellipsoid's volume.
                
                Buffer_frac(ind_list) = (1 - VolRatio);
                
                %BufferGrid(ind_list) = BufferGrid(ind_list)*(1 - VolRatio);
                %BufferVol = BufferVol + (1-VolRatio)*squeeze(sum(sum(ind_list~=0,1),2));
                
            end
            
            Elec_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount))...
                = Elec_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)) + Elec_frac;
                
            Buffer_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount))...
                = Buffer_grid(x_bot_indices(atomcount):x_top_indices(atomcount),y_bot_indices(atomcount):y_top_indices(atomcount), z_bot_indices(atomcount):z_top_indices(atomcount)).*Buffer_frac;
            
            clear Elec_frac;
            %toc
        end
        ind_list = Buffer_grid < 0;
        Buffer_grid(ind_list) = 0;
%         ind_list = BufferGrid == 1;
%         BufferGrid(ind_list) = 0;
        BufferVol = squeeze(sum(sum(Buffer_grid,1),2));
        
        toc
        
        ED_prof = zeros(zfgridnum,3);
        parray = zeros(zfgridnum,1);
        barray = zeros(zfgridnum,1);
        earray = zeros(zfgridnum,1);
        tote = 0;
        
        numgridarea = any(Elec_grid,3);
        minareagrid = sum(numgridarea(:));
        numemptygrid = sum(~numgridarea(:));
        
        for lay = 1:zfgridnum
                
            ED_prof(lay,1) = zf_grid(1,1,lay);
            %ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*xfgridnum*yfgridnum);
            ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*minareagrid);
            ED_prof(lay,3) = (BufferVol(lay)-numemptygrid)*fvol/(fvol*minareagrid);
            %ED_prof(lay,3) = BufferVol(lay)*fvol/(fvol*xfgridnum*yfgridnum);
            %ED_prof(lay,3) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay))+BufferVol(lay))*fvol/(fvol*xfgridnum*yfgridnum);

            parray(lay) = nnz(Elec_grid(:,:,lay))*fvol;
            barray(lay) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol;
            earray(lay) = sum(sum(Elec_grid(:,:,lay)));
            tote = tote + earray(lay);
        end
        
        ED_prof = flipud(ED_prof);
        for lay = 1:zfgridnum
            fprintf(wrf1,'%.8f %.8f %.8f\n', ED_prof(lay,1), ED_prof(lay,2), ED_prof(lay,3));
        end
        fclose('all');
        
        sprintf('minimal xy area of the box = %f',minareagrid*xfslice*yfslice)
        
        R.etotpdb = etotpdb;
        display(etotpdb);
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xyarea = abs(xtop-xbot)*abs(ytop-ybot);
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).minarea = minareagrid*xfslice*yfslice;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xyarea_protonly = area_nocov;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).etot = earray;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).tote = tote;
        display(sum(tote));
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).elecdif = (tote - etotpdb);
        display(sum(tote) - etotpdb);
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xfgridnum = xfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).yfgridnum = yfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).zfgridnum = zfgridnum;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xfgridnum_nocov = xfgridnum_nocov;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).yfgridnum_nocov = yfgridnum_nocov;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).zfgridnum_nocov = zfgridnum_nocov;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).fvol = fvol;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).fvol_nocov = fvol_nocov;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).provol = parray;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).bufvol = barray;
        R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).bufprot = BufferVol*fvol;
    end
end

end