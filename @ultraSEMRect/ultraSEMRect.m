classdef ultraSEMRect < ultraSEMQuad
%ULTRASEMRECT   Rectangular mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = ultraSEMRect( v )
            %ULTRASEMQUAD  Class constructor for the @ultraSEMQuad class.

            % Construct empty ultraSEMQuad:
            if ( nargin == 0 )
                return
            end
            
            if ( iscell(v) )
                for k = 1:numel(v)
                    obj(k,1) = ultraSEMRect( v{k} );
                end
                return
            end
            
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,2) ~= 2 ), v = v.'; end
            
            if ( size(v,2) == 2  && size(v, 1) == 4 )
                obj.v = v;
            else
                error('Incorrect vertices dimension.')
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function out = isRect(M)
            out = true;
        end
        
        function [T, mergeIdx] = refine(T, m)
            
            if ( nargin == 1 )
                m = 1;
            elseif ( m == 0 )
                mergeIdx = [];
                return
            end

            v = quad2rect(T.v);
            mergeIdx = {};
            
            for l = 1:m % Refine m times.  
                nDom = size(v, 1);
                dom = mat2cell(v, ones(nDom, 1)); % Convert to cell.
                for k = 1:nDom
                    d = dom{k}; % Subdivide this square into four new pieces:
                    midPt  = [ mean(d(1:2)), mean(d(3:4)) ];
                    dom{k} = [ d(1), midPt(1), d(3), midPt(2) ;
                        midPt(1), d(2), d(3), midPt(2) ;
                        midPt(1), d(2), midPt(2), d(4) ;
                        d(1), midPt(1), midPt(2), d(4)];
                end
                v = (vertcat(dom{:}));
                
                hMerge = reshape(1:4*nDom, 2, 2*nDom).';   % New horizontal merge.
                vMerge = reshape(1:2*nDom, 2, nDom).';     % New vertical merge.
                mergeIdx = [hMerge, vMerge, mergeIdx];     % Append to existing.
            end
            
            T = ultraSEMRect(rect2quad(v));
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end
