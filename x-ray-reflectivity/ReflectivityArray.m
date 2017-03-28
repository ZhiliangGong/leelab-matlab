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
        
    end
    
end