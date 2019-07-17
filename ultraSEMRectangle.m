classdef ultraSEMRectangle < ultraSEMQuad

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties ( Access = public )

        domain

    end

    
    methods
        
        function obj = ultraSEMRectangle( v )
            
            if ( nargin == 0 || isempty (v) )
                return
            end
            
            if ( size(v, 2) == 4 )
                n = size(v, 1);
                if ( n > 1 )
                    obj(n,1) = ultraSEMRectangle();
                    for k = 1:n
                        obj(k) = ultraSEMRectangle(v(k,:));
                    end
                else
                    obj.domain = v;
                end
                return
            end
            
%             obj.T1      = @(x,y) x;
%             obj.T2      = @(x,y) y;
%             obj.invT1   = @(x,y) x;
%             obj.invT2   = @(x,y) y;
%             obj.d2invT1 = @(x,y) 1+0*x;
%             obj.d2invT2 = @(x,y) 1+0*y;


            dom = [min(v(:,1)),max(v(:,1)), min(v(:,2)), max(v(:,2))];
            obj.domain = dom;
            
        end
        
        function v = vertices( r )
            r.domain
            v = r.domain([1 3 ; 2 3 ; 2 4 ; 1 4]);
        end
        
        function c = centroid( r )
            c = mean(r.domain([1 2 ; 3 4]), 2);
        end
        
    end
            
            
    
end


function v = four2eight(v)
    v = v([1 3 ; 2 3 ; 2 4 ; 1 4]);
end
    
    