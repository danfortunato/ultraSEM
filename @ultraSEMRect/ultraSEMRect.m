classdef ultraSEMRect < ultraSEMQuad
%ULTRASEMRECT   Rectangular mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = ultraSEMRect( v )
            %ULTRASEMRECT  Class constructor for the @ultraSEMRect class.

            % Construct empty ultraSEMRect:
            if ( nargin == 0 )
                return
            end
            
            % Construct multiple ultraSEMRects:
            if ( iscell(v) )
                obj(numel(v),1) = ultraSEMRect();
                for k = 1:numel(v)
                    obj(k,1) = ultraSEMRect( v{k} );
                end
                return
            end
            
            % Validate vertices:
%             v = ultraSEMRect.assertIsRect(v); % Necessary?
            
            % Assign vertices:
            obj.v = v;
            
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function out = isRect(R)
            out = true;
        end
        
        function v = rectVertices(R)
            v = quad2rect(R.v);
        end
        
        function v = quadVertices(R)
            v = rect2quad(R.v);
        end
        
        function [X, Y, XY] = transformGrid(T, x, y)
            rect = rectVertices(T);
            domx = rect(1:2);    domy = rect(3:4); 
            sclx = 2/diff(domx); scly = 2/diff(domy);
            X = (x+1)/sclx + domx(1); 
            Y = (y+1)/scly + domy(1);
            if ( nargout == 3 ), XY = [X Y]; end
        end
        
        function normal_d = transformNormalD(T, xy, p)
        %TRANSFORMNORMALD   Normal derivative operator for mapped domains.
        
            persistent bcrows_d
            if ( size(bcrows_d, 2) ~= p )
                % Construct normal derivatives conditions along the four edges:
                I = speye(p);
                lbc_d = kron( (-1).^(0:p-1).*(0:p-1).^2, I );
                rbc_d = kron( ones(1,p).*(0:p-1).^2, I );
                dbc_d = kron( I, (-1).^(0:p-1).*(0:p-1).^2 );
                ubc_d = kron( I, ones(1,p).*(0:p-1).^2 );
                bcrows_d = [ lbc_d ; rbc_d ; dbc_d ; ubc_d ];
            end

            rect = rectVertices(T);
            domx = rect(1:2);    domy = rect(3:4); 
            sclx = 2/diff(domx); scly = 2/diff(domy);
            normal_d = bcrows_d;
            normal_d(1:2*p,:) = sclx*normal_d(1:2*p,:);
            normal_d(2*p+1:end,:) = scly*normal_d(2*p+1:end,:);
        
        end
        
        function [op, rhs] = transformPDO(dom, op, rhs)

        end
        
        function [R, mergeIdx] = refine(R, m)
        %REFINE    Refine ultraSEMRect objects.    
            
            mergeIdx = {};
            % Parse inputs:
            if ( nargin == 1 )
                m = 1;
            elseif ( m == 0 )
                return
            end

            v = quad2rect(R.v);         % Convert quad representation to rect.
            
            for l = 1:m                 % Refine m times.  
                n = size(v, 1);         % Number of Rects.
                vNew = zeros(4*n, 4);   % Initialise new vertices.
                
                % Subdivide each rectangle into four new pieces and append:
                for k = 1:n      
                    vk = v(k,:); 
                    mid  = [ mean(vk(1:2)), mean(vk(3:4)) ];
                    vNew((k-1)*4+(1:4),:) = [ vk(1), mid(1), vk(3), mid(2) ;
                                              mid(1), vk(2), vk(3), mid(2) ;
                                              mid(1), vk(2), mid(2), vk(4) ;
                                              vk(1), mid(1), mid(2), vk(4) ];
                end
                v = vNew;
                hIdx = reshape(1:4*n, 2, 2*n).';   % New horizontal merge.
                vIdx = reshape(1:2*n, 2, n).';     % New vertical merge.
                mergeIdx = [hIdx, vIdx, mergeIdx]; %#ok<AGROW> Append to existing.  
            end
            
            R = ultraSEMRect(rect2quad(v));        % Assign new vertices..
            
        end
        
                    
%         function T = refinex(T, m)
%         %REFINEX   Refine a domain in the x-direction.
%         %   REFINEX(T) will divide each subdomain of T horizontally into
%         %   two new equally-sized pieces. The tree index information in the
%         %   result is updated to reflect the new subdomains, which are the
%         %   first to be merged.
%         
%         %   REFINEX(T, M) will refine M times.
% 
%             if ( isempty(T.domain) )
%                 return
%             end
%             if ( isempty(T.mergeIdx) && length(T) > 1)
%                 warning('Empty tree index encountered. Building a default one.');
%                 T.mergeIdx = defaultIdx(T.domain);
%             end
%             if ( nargin < 2 )
%                 m = 1;
%             end
% 
%             for l = 1:m % Refine m times.
%                 nDom = size(T.domain, 1);
%                 dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
%                 for k = 1:nDom
%                     d = dom{k}; % Subdivide this square into two new pieces:
%                     midPt  = mean(d(1:2));
%                     dom{k} = [ d(1), midPt, d(3:4) ;
%                                midPt, d(2), d(3:4) ];
%                 end
%                 T.domain = cell2mat(dom);              % Revert to a matrix.
%                 hMerge = reshape(1:2*nDom, 2, nDom).'; % New horizontal merge.
%                 T.mergeIdx = [hMerge, T.mergeIdx];     % Append to existing.
%             end
% 
%         end
% 
%         function T = refiney(T, m)
%         %REFINEY   Refine a domain in the y-direction.
%         %   REFINEY(T) will divide each subdomain of T vertically into two
%         %   new equally-sized pieces. The tree index information in the
%         %   result is updated to reflect the new subdomains, which are the
%         %   first to be merged.
%         %
%         %   REFINEY(T, M) will refine M times.
% 
%             if ( isempty(T.domain) )
%                 return
%             end
%             if ( isempty(T.mergeIdx) && length(T) > 1)
%                 warning('Empty tree index encountered. Building a default one.');
%                 T.mergeIdx = defaultIdx(T.domain);
%             end
%             if ( nargin < 2 )
%                 m = 1;
%             end
% 
%             for l = 1:m % Refine m times.
%                 nDom = size(T.domain, 1);
%                 dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
%                 for k = 1:nDom
%                     d = dom{k}; % Subdivide this square into two new pieces:
%                     midPt  = mean(d(3:4));
%                     dom{k} = [ d(1:2), d(3), midPt ;
%                                d(1:2), midPt, d(4) ];
%                 end
%                 T.domain = cell2mat(dom);              % Revert to a matrix.
%                 hMerge = reshape(1:2*nDom, 2, nDom).'; % New horizontal merge.
%                 T.mergeIdx = [hMerge, T.mergeIdx];     % Append to existing.
%             end
% 
%         end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = private, Static = true )
        
        function v = assertIsRect( v )
        %ASSERTISRECT    Check we have valid vertices for a ultraSEMRect.
            
            if ( ~isnumeric(v) )
                error('Input should be numeric.');
            end
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,2) ~= 2 ), v = v.'; end
            % Check dimension:
            if ( size(v,2) ~= 2 || size(v, 1) ~= 4 )
                error('Incorrect vertices dimension.')
            end
            % Check it's a valid rectangle:
            ndx = ~diff([v(:,1) ; v(1,1)]);
            ndy = ~diff([v(:,2) ; v(1,2)]);
            s1 = [0;1;0;1];
            s2 = [1;0;1;0];
            if ( ~( all(ndx==s1 & ndy==s2) || all(ndx==s2 & ndy==s1) ) )
                error('Not a valid rectangle')
            end
            
        end

    end

end
