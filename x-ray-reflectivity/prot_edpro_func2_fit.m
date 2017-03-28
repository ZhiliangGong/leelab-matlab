function [ED_prof] = prot_edpro_func2_fit(adata, theta, phi) %, currentFolder, pathname, pdbflag)
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

    %Now start to set the x,y,z axis

    %Rotate protein with theta and phi angles both between 0 to 2pi, where phi
    %is the rotation angle about z axis and theta is the rotation angle about x
    %axis


    % The phi rotation matrix about the Z axis following right hand rule is applied first using the martix:
    % [ cos(phi) -sin(phi) 0 ]
    % [ sin(phi)  cos(phi) 0 ]
    % [ 0         0        1 ]

    % The theta rotation matrix about the x axis following right hand rule is applied next using the matrix:
    % [ 1     0           0      ]
    % [ 0 cos(theta) -sin(theta) ]
    % [ 0 sin(theta)  costheta)  ]
    % 

    x = adata.x*cos(phi) - adata.y*sin(phi);
    y = adata.x*cos(theta)*sin(phi) + adata.y*cos(theta)*cos(phi) - adata.z*sin(theta);
    z = adata.x*sin(theta)*sin(phi) + adata.y*sin(theta)*cos(phi) + adata.z*cos(theta);


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
    %   3%%%%%%%%%%%%%%%%%%%4     %
        %%%%%%%%%%%%%%%%%%%%%     %
        %%%%%%%%%%%%%%%%%%%%%     %
        %%%%%%%%%%%%%%%%%%%%%     %
        %%%%%%%5%%%%%%%%%%%%%%%%%%6
        %%%%%%%%%%%%%%%%%%%%%%%%
    %   1%%%%%%%%%%%%%%%%%%%2

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

    for lay = 1:zfgridnum

        ED_prof(lay,1) = zf_grid(1,1,lay);
        ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*xfgridnum*yfgridnum);
        ED_prof(lay,3) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol/(fvol*xfgridnum*yfgridnum);

        parray(lay) = nnz(Elec_grid(:,:,lay))*fvol;
        barray(lay) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol;
        earray(lay) = sum(sum(Elec_grid(:,:,lay)));
        tote = tote + earray(lay);
    end

    ED_prof = flipud(ED_prof);

end