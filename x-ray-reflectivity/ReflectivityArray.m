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
        
        function overlayData(this, selected)
            
            if nargin == 1
                selected = 1 : length(this.data);
            end
            
            figure;
            hold on;
            
            legends = cell(1, length(selected));
            
            colorIndex = 0;
            for i = selected
                colorIndex = colorIndex + 1;
                this.plotData(this.data(i).proteinFit.refnorm, this.color(colorIndex));
                legends{colorIndex} = this.data(i).comment;
            end
            
            colorIndex = 0;
            for i = selected
                colorIndex = colorIndex + 1;
                [m, n] = this.bestFitIndices(this.data(i));
                this.plotFit(this.data(i).proteinFit.ref_fit{m, n}, this.color(colorIndex));
            end
            
            hold off;
            
            this.addAxesTitles();
            legend(legends);
            
        end
        
    end
    
    methods(Static)
        
        function plotDataAndFit(ref, refFit, color)
            
            if nargin == 2
                color = 'k';
            end
            
            errorbar(ref(:, 1), ref(:, 2), ref(:, 3), '.', 'color', color, 'linewidth', 1.2);
            plot(refFit(:, 1), refFit(:, 3), '-', 'color', color, 'linewidth', 1);
            
        end
        
        function plotData(ref, color)
            
            if nargin == 1
                color = 'k';
            end
            
            errorbar(ref(:, 1), ref(:, 2), ref(:, 3), '.', 'color', color, 'linewidth', 1.2);
            
        end
        
        function plotFit(ref, color)
            
            if nargin == 1
                color = 'k';
            end
            plot(ref(:, 1), ref(:, 3), '-', 'color', color, 'linewidth', 1);
            
        end
        
        function [m, n] = bestFitIndices(ref)
            
            index = find(ref.proteinFit.chi == min(ref.proteinFit.chi(:)), 1);
            [m, n] = ind2sub(size(ref.ed.profiles), index);
            
        end
        
        function color = color(n)
            
            colors = 'kbgcmry';
            index = mod(n, length(colors));
            if index == 0
                index = length(colors);
            end
            color = colors(index);
            
        end
        
        function addAxesTitles()
            
            set(gca, 'fontsize', 14);
            xlabel('$$ Q_z (\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('Normalized Reflectivity', 'fontsize', 16);
            
        end
        
    end
    
end