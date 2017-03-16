classdef DoubleLayer < handle
    
    properties
        
        surf
        thickness
        bulk
        decay
        
    end
    
    methods
        
        function this = DoubleLayer(surf, thickness, bulk, decay)
            
            this.surf = surf;
            this.thickness = thickness;
            this.bulk = bulk;
            this.decay = decay;
            
        end
        
        function distance = concPosition(this, conc)
            
            if conc > this.surf || conc < this.bulk
                distance = -1;
                warning('concentration must be equal to or smaller than the surf and bigger than or equal to the bulk');
            elseif conc == this.bulk
                distance = infinity;
            else
                distance = -log((conc - this.bulk) / this.surf) * this.decay;
            end
            
        end
        
        function value = surfAccumulation(this)
            
            value = this.surf * this.thickness * this.NA * 1e-28;
            
        end
       
        function value = bulkAccumulationWithinDistance(this, distance)
            
            value = integral(@(z) this.surf * exp(-z / this.decay) + this.bulk, 0, distance) * this.const;
            
        end
        
        function value = bulkAccumulationAtConc(this, conc)
            
            position = this.concPosition(conc);
            
            value = integral(@(z) this.surf * exp(-z / this.decay) + this.bulk, 0, position) * this.const;
            
        end
        
        function value = accumulationWithinDistance(this, distance)
            
            value = this.surfAccumulation + this.bulkAccumulationWithinDistance(distance);
            
        end
        
        function value = accumulationAtConc(this, conc)
            
            value = this.surfAccumulation + this.bulkAccumulationAtConc(conc);
            
        end
        
    end
    
    methods(Static)
        
        function NA = NA
            
            NA = 6.0221409e23;
            
        end
        
        function const = const
            
            const = 6.0221409e23 * 1e-28;
            
        end
        
    end
    
end