classdef TroughData
    %class for analyzing trough data
    %Zhiliang, Tiffany, Aug 31, 2016
    
    properties
        date
        lipid
        mw
        concentration
        volume
        subphase
        
        time
        areaPerMol
        surfPressure
        subTemp
        upTemp
        troughArea
        velocity
 
        operator
    end
    
    methods
        
        function x = TroughData(fname)
            
            if nargin == 0
                fname = uigetfile('*','Select the FM trough data to import.');
            end
            
            fid = fopen(fname);
            
            if ~fid
                error('Cannot open file.');
            end
            
            header = cell(10,1);
            for i = 1:10
                header{i} = fgetl(fid);
            end
            
            formatSpec = '%f %f %f %f %f %f %f';
            matrix = textscan(fid,formatSpec);
            
            fclose(fid);
            
            %read the header
            C = textscan(header{2},'%s %s %{MM/dd/uuuu}D %f %f %f %s','delimiter','\t');
            x.lipid = C{1}{1};
            x.operator = C{2}{1};
            x.date = C{3};
            x.volume = C{4};
            x.mw = C{5};
            x.concentration = C{6};
            x.subphase = C{7}{1};
            
            %process the matrix
            x.time = matrix{1};
            x.areaPerMol = matrix{2};
            x.surfPressure = matrix{3};
            x.subTemp = matrix{4};
            x.upTemp = matrix{5};
            x.troughArea = matrix{6};
            x.velocity = matrix{7};
            
        end
        
        function plotIsotherm(x,whole,formmatting)
            
            if whole
                n = length(x.time);
            else
                n = find(x.surfPressure==max(x.surfPressure),1);
            end
            
            if nargin == 2
                plot(x.areaPerMol(1:n),x.surfPressure(1:n),'k-','linewidth',2.4);
                xlabel('Area per Molecule (A^2)','fontsize',16);
                ylabel('\Pi (mN/m)','fontsize',16);
                set(gca,'fontsize',14);
            else
                if formmatting == 1
                    plot(x.areaPerMol,x.surfPressure,'k-','linewidth',2.4);
                    xlabel('Area per Molecule (A^2)','fontsize',16);
                    ylabel('\Pi (mN/m)','fontsize',16);
                    set(gca,'fontsize',14);
                else
                    plot(x.areaPerMol,x.surfPressure);
                end
            end
            
            
        end      
        
    end
    
end

