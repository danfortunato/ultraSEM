classdef Domain
%ULTRASEM.DOMAIN   Domain object from the ULTRASEM system.
%   An ULTRASEM.DOMAIN defines the domain on which we will solve a PDE.
%   This domain is formed (at least for now) by the union of the
%   rectangular domains stored in the .domain property. The ULTRASEM.DOMAIN
%   class is also responsible for merging patches during the build phase.
%   The merge information is contained in the .mergeIdx property (which is
%   documented in more detail below).

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain = []   % Domain of the patch.
        mergeIdx = {} % Determine which patches are merged at each level.
        % Each cell of mergeIdx denotes a level. Each row denotes the
        % indicies of two patches which are merged at that level. Each
        % patch index should occur exactly once per level. NaN means no
        % merge. The patch indices for the second level are determined by
        % the order of the first level, and so on.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Domain(dom, mergeIdx)

            % Construct empty domain:
            if ( nargin == 0 )
                return
            end

            if ( ischar(dom) )
                obj = obj.alphabet(dom);
                return
            end

            % Assign the domain:
            obj.domain = dom;

            % Assign the merge index (if given):
            if ( nargin > 1 )
                obj.mergeIdx = mergeIdx;
            else
                % Just make the default guess at the tree:
                obj.mergeIdx = ultraSEM.Domain.defaultIdx(dom);
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        % Construct alphabet domains:
        T = alphabet(str);

        % Make a default (naive) merge index for a domain.
        mergeIdx = defaultIdx(dom);

        % Construct rectangular domains:
        T = rectangle(varargin);

        % Construct quadrilateral domains:
        T = quad(varargin);

        % Construct triangular domains:
        T = triangle(varargin);

        % Construct triangular domains:
        T = duffy(varargin);

        % Construct convex polygonal domains:
        T = polygon(varargin);

        % Construct triangulated meshes:
        T = trimesh(varargin);

        % Construct word-shaped domains:
        T = scribble( str, join );

        % Check if a polygon is clockwise orientated:
        out = isClockwise( vertices );
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         function Q = refineCornerRect(d, k)
% %             d = T.domain;
%             c = mean(d([1 3 ; 2 4]))';
%             v = d([1 2 2 1 ; 3 3 4 4]);
%             if ( k == 1 )
%                 vnew = (v(:,1) + v(:,[2,4]))/2;
%                 v1 = [v(:,1) vnew(:,1) c vnew(:,2)];
%                 v2 = [v(:,[2 3]) c vnew(:,1)];
%                 v3 = [v(:,4) vnew(:,2) c v(:,3)];
%             elseif ( k == 2 )
%                 vnew = (v(:,2) + v(:,[3,1]))/2;
%                 v1 = [v(:,2) vnew(:,1) c vnew(:,2)];
%                 v2 = [v(:,[3 4]) c vnew(:,1)];
%                 v3 = [v(:,1) vnew(:,2) c v(:,4)];                
%             elseif ( k == 3 )
%                 vnew = (v(:,3) + v(:,[4,2]))/2;
%                 v1 = [v(:,3) vnew(:,1) c vnew(:,2)];
%                 v2 = [v(:,[4 1]) c vnew(:,1)];
%                 v3 = [v(:,2) vnew(:,2) c v(:,1)]; 
%             elseif ( k == 4 )
%                 vnew = (v(:,4) + v(:,[1,3]))/2;
%                 v1 = [v(:,4) vnew(:,1) c vnew(:,2)];
%                 v2 = [v(:,[1 2]) c vnew(:,1)];
%                 v3 = [v(:,3) vnew(:,2) c v(:,2)];                 
%             end            
%             Q(3,1) = ultraSEM.Quad();
%             Q(1) = ultraSEM.Quad(v1);
%             Q(2) = ultraSEM.Quad(v2);
%             Q(3) = ultraSEM.Quad(v3);
%         end

