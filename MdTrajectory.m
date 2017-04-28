classdef MdTrajectory < handle
    
    properties
        
        traj
        mass
        time
        
        best_frame
        rmsd_mean
        
    end
    
    methods
        
        function this = MdTrajectory(file, time)
            
            this.traj = readdcd(file);
            
            if nargin == 1
                time = 1 : size(this.traj, 1);
            end
            
            this.time = time;
            
        end
        
        function getMeanRmsd(this)
            
            n = this.length();
            
            trj = this.traj;
            mat = zeros(n, n);
            parfor i = 1 : n
                mat(:, i) = superimpose(trj(i, :), trj);
            end
            
            this.rmsd_mean = mean(mat, 2)';
            
        end
        
        function n = getBestFrame(this)
            
            if isempty(this.rmsd_mean)
                this.getMeanRmsd();
            end
            
            n = find(this.rmsd_mean == min(this.rmsd_mean));
            
        end
        
        function n = getMostSimilarFrame(this, trj)
            
            rmsd = superimpose(trj, this.traj);
            n = find(rmsd == min(rmsd), 1);
            
        end
        
        % plot
        
        function plotMeanRmsd(this)
            
            if isempty(this.rmsd_mean)
                this.getMeanRmsd();
            end
            
            figure;
            plot(this.time, this.rmsd_mean, '-', 'linewidth', 2.4);
            xlabel('Time (ns)', 'fontsize', 16);
            ylabel('RMSD (\AA)', 'fontsize', 16, 'interpreter', 'latex');
            set(gca, 'fontsize', 14);
            
        end
        
        % utility
        
        function n = length(this)
            
            n = size(this.traj, 1);
            
        end
        
    end
    
    methods (Static)
        
        
    end
    
end