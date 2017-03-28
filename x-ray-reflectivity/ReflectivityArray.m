classdef ReflectivityArray < handle
    
    properties
        
        data
        
    end
    
    methods
        
        function this = ReflectivityArray(ref)
            
            if nargin == 1
               this.push(ref); 
            end
            
        end
        
        function push(this, ref)
            
            if isa(ref, 'Reflectivity')
                if isempty(this.data)
                    this.data = ref;
                else
                    this.data(end + 1) = ref;
                end
            end
            
        end
        
        function overlayData(this)
            
            figure;
            hold on;
            for i = 1 : length(this.data)
                [m, n] = this.bestFitIndices(this.data(i));
                this.plotDataAndFit(this.data(i).proteinFit.refnorm, this.data(i).proteinFit.ref_fit{m, n});
            end
            hold off;
            
        end
        
    end
    
    methods(Static)
        
        function plotDataAndFit(ref, refFit)
            
            errorbar(ref(:, 1), ref(:, 2), ref(:, 3), '.', 'color', [0, 0, 0] + 0.5, 'linewidth', 1.2);
            plot(refFit(:, 1), refFit(:, 3), '-k', 'linewidth', 2.4);
            set(gca, 'fontsize', 14);
            xlabel('$$ Q_z (\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('Normalized Reflectivity', 'fontsize', 16);
            legend('Data', 'Fit');
            
        end
        
        function [m, n] = bestFitIndices(ref)
            
            index = find(ref.proteinFit.chi == min(ref.proteinFit.chi(:)), 1);
            [m, n] = ind2sub(size(ref.ed.profiles), index);
            
        end
        
    end
    
end