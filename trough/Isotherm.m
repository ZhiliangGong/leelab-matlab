classdef Isotherm < handle
    
    properties
        raw
        time
        area
        pressure
        molArea
    end
    
    methods
        
        function this = Isotherm(file)
            
            if nargin == 0
                [file, path] = uigetfile('*', 'Select an isotherm file');
                file = fullfile(path, file);
            end
            
            try
                this.raw = importdata(file);
                this.time = this.raw.data(:, 1);
                this.area = this.raw.data(:, 2);
                this.molArea = this.raw.data(:, 3);
                this.pressure = this.raw.data(:, 6);
            catch EM
                disp(EM);
            end
            
        end
        
        function plotTimeArea(this, varargin)
            
            startTime = min(this.time);
            finishTime = max(this.time);
            
            if ~isempty(varargin) > 0
                startTime = varargin{1};
                if length(varargin) > 1
                    finishTime = varargin{2};
                end
            end
            
            startIndex = find((this.time - startTime).^2 == min((this.time - startTime).^2), 1);
            finishIndex = find((this.time - finishTime).^2 == min((this.time - finishTime).^2), 1);
            
            figure;
            plot(this.time(startIndex : finishIndex), this.area(startIndex : finishIndex), 'k-', 'linewidth', 2.4);
            set(gca, 'fontsize', 14);
            xlabel('Time (s)', 'fontsize', 16);
            ylabel('Area (cm^2)', 'fontsize', 16);
        end
        
        function plotIsotherm(this, varargin)
            
            startIndex = 1;
            finishIndex = length(this.time);
            
            if ~isempty(varargin) > 0
                startTime = varargin{1};
                startIndex = find((this.time - startTime).^2 == min((this.time - startTime).^2), 1);
                if length(varargin) > 1
                    finishTime = varargin{2};
                    finishIndex = find((this.time - finishTime).^2 == min((this.time - finishTime).^2), 1);
                end
            end
            
            figure;
            plot(this.molArea(startIndex : finishIndex), this.pressure(startIndex : finishIndex), 'k-', 'linewidth', 2.4);
            set(gca, 'fontsize', 14);
            xlabel('$$ Area / Molecule (\AA^2) $$', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('Surface Pressure (mN/m)', 'interpreter', 'latex', 'fontsize', 16);
            
        end
        
    end
    
end