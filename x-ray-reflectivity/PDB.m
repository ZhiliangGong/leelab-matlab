classdef PDB < handle
    
    properties
        
        file
        data
        ed
        
    end
    
    methods
        
        function this = PDB(file)
            
            if nargin == 0
                [filename, pathname] = uigetfile('*', 'Select the pdb file.');
                file = fullfile(pathname, filename);
            end
            
            this.file = file;
            text = fileread(file);
            this.data = processPdbFile(text);
            
            % get the mass and electrons of each atom as a vector
            n = length(this.data.x);
            this.data.mass = zeros(1, n);
            this.data.electron = zeros(1, n);
            this.data.radius = zeros(1, n);
            [massTable, electronTable] = this.getPeriodicTable();
            radiusTable = this.getAtomRadiusTable();
            for i = 1 : n
                this.data.mass(i) = massTable(this.data.atoms(i));
                this.data.electron(i) = electronTable(this.data.atoms(i));
                this.data.radius(i) = radiusTable(this.data.atoms(i));
            end
            
        end
        
        function generateEdProfiles(this, theta, phi)
            
            if nargin == 1
                theta = 0 : 5 : 180;
                phi = 0 : 10 : 350;
            end
            
            gridSize = 0.5;
            this.ed.gridSize = gridSize;
            this.ed.theta = theta;
            this.ed.phi = phi;
            
            m = length(theta);
            n = length(phi);
            this.ed.profiles = cell(m, n);
            
            for i = 1 : m
                for j = 1 : n
                    
                    positions = [1, 0, 0; 0, cos(theta(i)), -sin(theta(i)); 0, sin(theta(i)), cos(theta(i))]...
                        * [cos(phi(j)), -sin(phi(j)), 0; sin(phi(j)), cos(phi(j)), 0; 0, 0, 1]...
                        * [this.data.x; this.data.y; this.data.z];
                    x = positions(1,:);
                    y = positions(2,:);
                    z = positions(3,:);
                    
                    xtop = max(x + this.data.radius) + gridSize;
                    xbot = min(x - this.data.radius) - gridSize;
                    ytop = max(y + this.data.radius) + gridSize;
                    ybot = min(y - this.data.radius) - gridSize;
                    ztop = max(z + this.data.radius) + gridSize;
                    zbot = min(z - this.data.radius) - gridSize;
                    
                    
                    xfgridnum = round((xtop - xbot) / gridSize);
                    yfgridnum = round((ytop - ybot) / gridSize);
                    zfgridnum = round((ztop - zbot) / gridSize);
                    
                    xfslice = (xtop - xbot) / xfgridnum;
                    yfslice = (ytop - ybot) / yfgridnum;
                    zfslice = (ztop - zbot) / zfgridnum;
                    fvol = xfslice * yfslice * zfslice;
                    
                    Elec_grid = zeros(xfgridnum, yfgridnum, zfgridnum);
                    
                    x_grid_pos = xbot + xfslice * ((1:xfgridnum) - 0.5);
                    y_grid_pos = ybot + yfslice * ((1:yfgridnum) - 0.5);
                    z_grid_pos = zbot + zfslice * ((1:zfgridnum) - 0.5);
                    
                    x_diff = repmat(x_grid_pos, length(x), 1) - repmat(x', 1, length(x_grid_pos));
                    y_diff = repmat(y_grid_pos, length(y), 1) - repmat(y', 1, length(y_grid_pos));
                    z_diff = repmat(z_grid_pos, length(z), 1) - repmat(z', 1, length(z_grid_pos));
                    
                    [x_top_atom_indices(:,1), x_top_atom_indices(:,2)] = ind2sub([length(x), length(x_grid_pos)], find(x_diff > this.data.radius' & x_diff < this.data.radius' + xfslice));
                    [x_bot_atom_indices(:,1), x_bot_atom_indices(:,2)] = ind2sub([length(x), length(x_grid_pos)], find(x_diff < -this.data.radius' & x_diff > -this.data.radius' - xfslice));
                    [y_top_atom_indices(:,1), y_top_atom_indices(:,2)] = ind2sub([length(y), length(y_grid_pos)], find(y_diff > this.data.radius' & y_diff < this.data.radius' + yfslice));
                    [y_bot_atom_indices(:,1), y_bot_atom_indices(:,2)] = ind2sub([length(y), length(y_grid_pos)], find(y_diff < -this.data.radius' & y_diff > -this.data.radius' - yfslice));
                    [z_top_atom_indices(:,1), z_top_atom_indices(:,2)] = ind2sub([length(z), length(z_grid_pos)], find(z_diff > this.data.radius' & z_diff < this.data.radius' + zfslice));
                    [z_bot_atom_indices(:,1), z_bot_atom_indices(:,2)] = ind2sub([length(z), length(z_grid_pos)], find(z_diff < -this.data.radius' & z_diff > -this.data.radius' - zfslice));
                    
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
                    
                    tic
                    for k = 1:length(this.data.atoms)
                        
                        Elec_frac = zeros(x_top_indices(k)-x_bot_indices(k)+1, y_top_indices(k)-y_bot_indices(k)+1, z_top_indices(k)-z_bot_indices(k)+1);
                        
                        atomdist = (xf_grid(x_bot_indices(k):x_top_indices(k),y_bot_indices(k):y_top_indices(k), z_bot_indices(k):z_top_indices(k)) - x(k)).^2 ...
                            + (yf_grid(x_bot_indices(k):x_top_indices(k),y_bot_indices(k):y_top_indices(k), z_bot_indices(k):z_top_indices(k)) - y(k)).^2 ...
                            + (zf_grid(x_bot_indices(k):x_top_indices(k),y_bot_indices(k):y_top_indices(k), z_bot_indices(k):z_top_indices(k)) - z(k)).^2;
                        
                        ind_list = (atomdist <= this.data.radius(k)^2);
                        
                        Elec_frac(ind_list) = (this.data.electron(k)*fvol)/((4/3)*pi*this.data.radius(k)^3);
                        
                        Elec_grid(x_bot_indices(k):x_top_indices(k),y_bot_indices(k):y_top_indices(k), z_bot_indices(k):z_top_indices(k))...
                            = Elec_grid(x_bot_indices(k):x_top_indices(k),y_bot_indices(k):y_top_indices(k), z_bot_indices(k):z_top_indices(k)) + Elec_frac;
                        
                    end
                    
                    toc
                    
                    ED_prof = zeros(zfgridnum,3);
                    parray = zeros(zfgridnum,1);
                    barray = zeros(zfgridnum,1);
                    earray = zeros(zfgridnum,1);
                    tote = 0;
                    
                    numgridarea = any(Elec_grid,3);
                    minareagrid = sum(numgridarea(:));
                    
                    for lay = 1:zfgridnum
                        
                        ED_prof(lay,1) = zf_grid(1,1,lay);
                        ED_prof(lay,2) = sum(sum(Elec_grid(:,:,lay)))/(fvol*xfgridnum*yfgridnum);
                        ED_prof(lay,3) = sum(sum(Elec_grid(:,:,lay)))/(fvol*minareagrid);
                        ED_prof(lay,4) = (xfgridnum*yfgridnum-nnz(Elec_grid(:,:,lay)))*fvol/(fvol*xfgridnum*yfgridnum);
                        ED_prof(lay,5) = (minareagrid-nnz(Elec_grid(:,:,lay)))*fvol/(fvol*minareagrid);
                        
                        parray(lay) = nnz(Elec_grid(:,:,lay))*fvol;
                        barray(lay) = (minareagrid-nnz(Elec_grid(:,:,lay)))*fvol;
                        earray(lay) = sum(sum(Elec_grid(:,:,lay)));
                        tote = tote + earray(lay);
                    end
                    
                    ED_profout.(sprintf('t%03dp%03d',theta_rot,phi_rot)) = flipud(ED_prof);
                    R.etotpdb = etotpdb;
                    display(etotpdb);
                    R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).xyarea = abs(xtop-xbot)*abs(ytop-ybot);
                    R.(sprintf('t%03dp%03d',theta_rot,phi_rot)).minarea = minareagrid*xfslice*yfslice;
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
        
        % utility
        
        function total = totalElectron(this)
            
            total = sum(this.data.electron);
            
        end
        
        function mw = molecularWeight(this)
            
            mw = sum(this.data.mass);
            
        end
        
    end
    
    methods(Static)
        
        function [AtomicMassTable, ElectronTable] = getPeriodicTable()

            symbols = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U';
            atomicNumber = (1:92);
            masses = [1.0079,4.0026,6.941,9.0122,10.81,12.011,14.007,15.999,18.998,20.18,22.99,24.305,26.982,28.085,30.974,32.066,35.453,39.948,39.098,40.078,44.956,47.867,50.941,51.996,54.938,55.845,58.933,58.693,63.546,65.39,69.723,72.61,74.922,78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,98,101.07,102.91,106.42,107.87,112.41,114.82,118.71,121.76,127.6,126.9,131.29,132.91,137.33,138.91,140.12,140.91,144.24,145,150.36,151.96,157.25,158.93,162.5,164.93,167.26,168.93,173.04,174.97,178.49,180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,204.38,207.2,208.98,209,210,222,223,226,227,232.04,231.04,238.03];

            symbol = regexp(symbols, ' ', 'split');
            ElectronTable = containers.Map;
            AtomicMassTable = containers.Map;

            for i = 1:length(symbol)
                ElectronTable(symbol{i}) = atomicNumber(i);
                AtomicMassTable(symbol{i}) = masses(i);
            end

        end
        
        function radius = getAtomRadiusTable()
            
            radius = containers.Map;
            radius('C') = 1.7;
            radius('N') = 1.55;
            radius('H') = 1.1;
            radius('O') = 1.52;
            radius('S') = 1.8;
            radius('P') = 1.95;
            
        end
        
    end
    
end