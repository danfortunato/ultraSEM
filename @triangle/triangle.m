classdef triangle < ultraSEMDomain
%TRIANGLE   Kite mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        corners

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function T = refine(T, r)
            
            if ( nargin == 1 )
                r = 1;
            end
       
            v = T.corners;
            vnew = (v(1:3,:) + v([2 3 1],:))/2;
            T = {};
            T{1} = ultraSEMDomain.triangle([v(1,:) ; vnew(1,:) ; vnew(3,:)]);
            T{2} = ultraSEMDomain.triangle(vnew);
            T{3} = ultraSEMDomain.triangle([vnew(1,:) ; v(2,:) ; vnew(2,:)]);
            T{4} = ultraSEMDomain.triangle([vnew(2,:) ; v(3,:) ; vnew(3,:)]);
            
            if ( r > 1 )  
                T = cellfun(@(T) refine(T,r-1), T, 'uniformoutput', false);
            end
            T = merge(T{:});
%             T = T{1} & T{2} & T{3} & T{4};
            
        end
        
    end
    
end