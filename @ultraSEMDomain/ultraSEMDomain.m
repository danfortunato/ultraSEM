classdef ultraSEMDomain
%ULTRASEMDOMAIN  Domain object from the ULTRASEM system.
%   An ULTRASEMDOMAIN defines the domain on which we will solve a PDE. This
%   domain is formed (at least for now) by the union of the rectangular
%   domains stored in the .domain property. The ULTRASEMDOMAIN class is
%   also responsible for merging patches during the build phase. The merge
%   information is contained in the .mergeIdx property (which is documented
%   in more detail below).

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>

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
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEMDomain(dom, tree)
            %ULTRASEMDOMAIN   Class constructor for the @ultraSEMDomain class.

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

            % Assign the tree (if given):
            if ( nargin > 1 )
                obj.mergeIdx = tree;
            else
                % Just make the default guess at the tree:
                obj.mergeIdx = ultraSEMDomain.defaultIdx( dom );
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function C = and(A, B)
        %AND   Merge two ultraSEMDomains.
        %   C = A & B will merge the two ultraSEMDomains A and B. See
        %   MERGE() for further details.
        %
        % See also MERGE().
            C = merge(A, B); % This method is simply a wrapper for MERGE().
        end

        function p = build(T, p)

            % TODO: Document.

            if ( length(T) == 1 || numel(p) == 1 )
                % Nothing to build for a single patch.
                return
            elseif ( isempty(T.mergeIdx) )
                warning('ULTRASEM:ULTRASEMDOMAIN:build:noMergeIdx', ...
                    '%s does not contain merge information.', inputname(1))
                T.mergeIdx = T.defaultIdx(T.domain);
            end

            % Loop over each level:
            for j = 1:numel(T.mergeIdx) % <-- number of levels

                idxj = T.mergeIdx{j};         % Indicies for jth level.
                q = cell(size(idxj, 1),1);    % Initialise patches at jth level.

                % Perform all the merges at this level:
                for k = 1:size(idxj, 1)
                    idxjk = idxj(k,:);        % Index for kth merge at jth level.
                    idxjk(isnan(idxjk)) = []; % Remove NaNs.
                    pk = p(idxjk);            % Extract patches.
                    pk(cellfun(@isempty, pk)) = []; % Remove empty patches.
                    if ( isempty(pk) )
                        q{k} = [];            % Empty merge.
                    else
                        q{k} = merge(pk{:});  % Perform merge.
                    end
                end

                p = q; % Store back in p for recursion.

            end

        end

        function T = ctranspose(T)
        %CTRANSPOSE   Rotate an ultraSEMDomain clockwise by 90 degrees.
        %
        % See also TRANSPOSE().
            T = transpose(T); % This method is simply a wrapper for TRANSPOSE().
        end

        function T = fliplr(T)
        %FLIPLR   Flip an ultraSEMDomain horizontally about its center.
            minx = min(T.domain(:,1));
            maxx = max(T.domain(:,2));
            T.domain(:,[1 2]) = minx + maxx - T.domain(:,[2 1]);
        end

        function T = flipud(T)
        %FLIPUD   Flip an ultraSEMDomain vertically about its center.
            miny = min(T.domain(:,3));
            maxy = max(T.domain(:,4));
            T.domain(:,[3 4]) = miny + maxy - T.domain(:,[4 3]);
        end

        function [d, ia, ib] = intersect(S, T)
        %INTERSECT   Compute the intersecting patches of two ultraSEMDomains.
            [d, ia, ib] = intersect(S.domain, T.domain, 'rows', 'stable');
        end

        function out = length(T)
        %LENGTH   The length of an ultraSEMDomain is its number of patches.
            out = size(T.domain, 1);
        end

        function H = merge(varargin)
        %MERGE   Merge two or more ultraSEMDomains.
        
            F = varargin{1};
            G = varargin{2};
            
            % Merge multiple pieces:
            if ( nargin > 2 )
                H = merge(merge(F,G), varargin{3:end});
                return
            end
            % Nothing to do here:
            if ( isempty(F) )
                H = G;
                return
            elseif ( isempty(G) )
                H = F;
                return
            end

            treeSizeF = size(F.mergeIdx, 2);
            treeSizeG = size(G.mergeIdx, 2);
            % The output will have one more level than F or G.
            treeSizeH = max(treeSizeF, treeSizeG) + 1;

            % Pad F or G so that they have the same number of levels:
            if ( treeSizeF < treeSizeG )
                F.mergeIdx(end+1:treeSizeG) = ...
                    repmat({[1, NaN]}, 1, treeSizeG-treeSizeF);
            elseif ( treeSizeF > treeSizeG )
                G.mergeIdx(end+1:treeSizeF) = ...
                    repmat({[1, NaN]}, 1, treeSizeF-treeSizeG);
            end

            % To update the mrge information we need to adjust the merge
            % indicies of G by the number of patches in F for each level.
            % skip(k) contains this number.
            skip = cellfun(@(t) max(t(:)), F.mergeIdx);
            mergeIdx = cell(1, treeSizeH); % Initialise.
            for k = 1:(treeSizeH-1)
                sF = size(F.mergeIdx{k},2);
                sG = size(G.mergeIdx{k},2);
                if ( sF > sG )
                    G.mergeIdx{k}(:,end+1) = NaN;
                elseif ( sF < sG )
                    F.mergeIdx{k}(:,end+1) = NaN;
                end
                mergeIdx{k} = [F.mergeIdx{k}; G.mergeIdx{k} + skip(k)];
            end
            mergeIdx{end} = [1, 2]; % Final merge leaves a single patch.

            if ( isnumeric(F.domain) && ~isnumeric(G.domain) )
                F.domain = ultraSEM.rectangle(F.domain, 'makeObj');
            elseif  ( ~isnumeric(F.domain) && isnumeric(G.domain) )
                G.domain = ultraSEM.rectangle(G.domain, 'makeObj');
            end

            % Construct the new ultraSEMDomain:
            H = ultraSEMDomain([F.domain ; G.domain], mergeIdx);

        end

        function T = minus(T, c)
        %MINUS   Shift an ultraSEMDomain or setminus two ultraSEMDomains.
        %   T - C, where T is a ultraSEMDomain and C is a scalar, is
        %   equivalent to T + (-C).
        %
        %   T - S, where S and T are both ultraSEMDomains, will remove
        %   common patches of S and T from T.
        %
        % See also PLUS.

            if ( isnumeric(T) )
                T = plus(c, (-T));
            elseif ( isnumeric(c) )
                T = plus(T, (-c));
            elseif ( isa(T,'ultraSEMDomain') && isa(c,'ultraSEMDomain') )
                % TODO: Allow us to cut out shapes from a domain.
                [~, iT, ~] = intersect(T, c);
                T = removePatch(T, iT);
            else
                error('ULTRASEM:ULTRASEMDOMAIN:plus:unknown', ...
                    'Cannot subtract a %s from a %s.', class(c), class(T));
            end

        end

        function T = mtimes(T, c)
        %MTIMES   Scale an ultraSEMDomain.
        %   C*T will shift the ultraSEMDomain to the right by real(c) and
        %   upwards by imag(C). C must be a scalar.

            if ( ~isa(T, 'ultraSEMDomain') )
                % Ensure T is the ultraSEMDomain:
                T = mtimes(c, T);
                return
            elseif ( isa(c, 'ultraSEMDomain' ) )
                % We can't multiply two ultraSEMDomains:
                error('ULTRASEM:ULTRASEMDOMAIN:mtimes:twodomains', ...
                    'Cannot multiply (* or .*) two ultraSEMDomains.\n')
            elseif ( ~isnumeric(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:mtimes:unknown', ...
                    'Cannot multiply an object of type %s by an ultraSEMDomain.', ...
                    class(c));
            elseif ( ~isscalar(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:mtimes:scalar', ...
                    'C must be a scalar.')
            end

            % Scale the domain:
            T.domain = c*T.domain;

        end

        function E = not(T, pad)
        %NOT   Take the complement of an ultraSEMDomain.
        %   ~T forms a rectangular domain by constructing a rectangular
        %   domain R containing T and removing T from it. The size of R
        %   will be one patch size greater than the extremities of T.
        %
        %   NOT(T, P) is similar to the above, but will pad by an amount P
        %   on each side of T. If P is a 4x1 vector then it is interpreted
        %   as P = [P_left, P_right, P_bottom, P_top].

            % Determine padding size:
            if ( nargin == 1 )
                pad = [1 1 1 1];
            elseif ( isscalar(pad) )
                pad = repmat(pad, 1, 4);
            end

            % Determine sizes:
            dom = T.domain;
            dx = diff(dom(1,1:2)); dy = diff(dom(1,3:4));
            x1 = min(dom(:,1));    x2 = max(dom(:,2));
            y1 = min(dom(:,3));    y2 = max(dom(:,4));

            % Construct containing rectangle:
            newDom = [x1-pad(1)*dx, x2+pad(2)*dx y1-pad(3)*dy y2+pad(4)*dy];
            m = diff(newDom(1:2))/dx; n = diff(newDom(3:4))/dy;
            R = ultraSEM.rectangle(newDom, m, n);

            % Extract T from R:
            E = R - T;

        end

        function varargout = plot(T, varargin)
        %PLOT   Plot an ultraSEMDomain.
        %   PLOT(T) plots the ultraSEMDomain T and indicates the patch
        %   numbering for the first layer.
        %
        %   PLOT(T, C), where C is a single character string chosen from
        %   the list 'r', 'g', 'b', 'c', 'm', 'y', 'w', 'k', or an RGB row
        %   vector triple, [r g b], fills the domain with the constant
        %   specified color.
        %
        %   PLOT(T, C, PROP1, VAL1, PROP2, VAL2) allows adjusting of the
        %   plot in any way allowed by the built-in FILL() method.
        %
        %   H = PLOT(T, ...) returns a figure handle of the form returned
        %   by H = FILL(...), where FILL() is the built-in MATLAB method.

            if ( ~isnumeric(T.domain) )
                [varargout{1:nargout}] = plot(T.domain, varargin{:});
                return
            end

            % Choose a color:
            if ( nargin > 1 && ischar(varargin{1}) && ...
                    ~isempty(regexp( varargin{1}, '[bgrcmykw]', 'match')) )
                c = varargin{1};
                varargin(1) = [];
            elseif ( nargin > 1 && isnumeric(varargin{1}) )
                c = varargin{1};
                varargin(1) = [];
            else
                c = 'm';
            end

            holdState = ishold();

            % Loop over the patches:
            for k = 1:length(T)
                domk = T.domain(k,:);
                % Call the built-in FILL method:
                h(k,1) = fill(domk([1 2 2 1]), domk([3 3 4 4]), c, ...
                    'FaceAlpha', .25, varargin{:}); hold on %#ok<AGROW>
                % Add text to center of patch:
                if ( length(T) < 100 ) % Don't add text on large meshes
                    text(mean(domk(1:2)), mean(domk(3:4)), int2str(k), ...
                        'HorizontalAlignment', 'center')
                end
            end
            axis equal

            % Return hold state:
            if ( ~holdState), hold off, end
            % Assign output if required:
            if ( nargout > 0 ), varargout = {h}; end

        end

        function T = plus(T, c)
        %PLUS   Shift an ultraSEMDomain.
        %   T + C will shift the ultraSEMDomain to the right by real(c) and
        %   upwards by imag(C). C must be a scalar.
        %
        % See also MINUS().

            if ( ~isa(T, 'ultraSEMDomain') )
                % Ensure T is the ultraSEMDomain:
                T = plus(c, T);
                return
            elseif ( isa(c, 'ultraSEMDomain' ) )
                % We can't add two ultraSEMDomains:
                error('ULTRASEM:ULTRASEMDOMAIN:plus:twodomains', ...
                    ['Cannot add (+) two ultraSEMDomains.\n', ...
                       'Did you mean to merge them with &?'])
            elseif ( ~isnumeric(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:plus:unknown', ...
                    'Cannot add an object of type %s to a ultraSEMDomain.', ...
                    class(c));
            elseif ( ~isscalar(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:plus:scalar', ...
                    'C must be a scalar.')
            end

            % Shift the domain:
            T.domain(:,1:2) = T.domain(:,1:2) + real(c);
            T.domain(:,3:4) = T.domain(:,3:4) + imag(c);

        end

        function T = refine(T, m)
        %REFINE   Refine an ultraSEMDomain.
        %   REFINE(T) will divide each subdomain of T into four new
        %   equally-sized pieces. The tree index information in the result
        %   is updated to reflect the new subdomains, which are the first
        %   to be merged.
        %
        %   REFINE(T, M) will refine M times.

            if ( isempty(T.domain) )
                return
            end
            if ( isempty(T.mergeIdx) && length(T) > 1)
                warning('Empty tree index encountered. Building a default one.');
                T.mergeIdx = defaultIdx(T.domain);
            end
            if ( nargin < 2 )
                m = 1;
            end
            
            if ( isnumeric(T.domain) )
                T = refineRectangle(T, m);
            else
                T= refineQuad(T,m);
            end
                

        end
        
        function T = refineRectangle(T, m)
            for l = 1:m % Refine m times.
                nDom = size(T.domain, 1);
                dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
                for k = 1:nDom
                    d = dom{k}; % Subdivide this square into four new pieces:
                    midPt  = [ mean(d(1:2)), mean(d(3:4)) ];
                    dom{k} = [ d(1), midPt(1), d(3), midPt(2) ;
                        midPt(1), d(2), d(3), midPt(2) ;
                        d(1), midPt(1), midPt(2), d(4) ;
                        midPt(1), d(2), midPt(2), d(4) ];
                end
                T.domain = cell2mat(dom);                  % Revert to a matrix.
                hMerge = reshape(1:4*nDom, 2, 2*nDom).';   % New horizontal merge.
                vMerge = reshape(1:2*nDom, 2, nDom).';     % New vertical merge.
                T.mergeIdx = [hMerge, vMerge, T.mergeIdx]; % Append to existing.
            end
        end
            
        
        function T2 = refineQuad(T, m)
            if ( m > 1 )
                T2 = refine(refine(T,1), m-1);
                return
            end
            
            nDom = size(T.domain, 1);
            T2 = cell(nDom, 1);
            for k = 1:nDom
                Qk = T.domain(k);
                xRef = Qk.domain([1 2 2 1]);
                yRef = Qk.domain([3 3 4 4]);
                x = Qk.T1(xRef, yRef);
                y = Qk.T2(xRef, yRef);
                v = [x; y];
                
                v2 = (v(:,[1 2 3 4]) + v(:,[2 3 4 1]))/2;
                
                p = polyshape(x, y);
                [xmid, ymid] = centroid(p);
                midpt = [xmid ; ymid];
                
                Qk2(4,1) = quad();
                Qk2(1) = quad([v(:,1) v2(:,1) midpt v2(:,4)]');
                Qk2(2) = quad([v(:,2) v2(:,2) midpt v2(:,1)]');
                Qk2(3) = quad([v(:,3) v2(:,3) midpt v2(:,2)]');
                Qk2(4) = quad([v(:,4) v2(:,4) midpt v2(:,3)]');
                
                T2{k} = ultraSEMDomain(Qk2, {[1 2 ; 3 4], [1 2]});
            end
            T2 = merge(T2{:});
        end
 
        

        function T = refinex(T, m)
        %REFINEX   Refine a domain in the x-direction.
        %   REFINEX(T) will divide each subdomain of T horizontally into
        %   two new equally-sized pieces. The tree index information in the
        %   result is updated to reflect the new subdomains, which are the
        %   first to be merged.
        
        %   REFINEX(T, M) will refine M times.

            if ( isempty(T.domain) )
                return
            end
            if ( isempty(T.mergeIdx) && length(T) > 1)
                warning('Empty tree index encountered. Building a default one.');
                T.mergeIdx = defaultIdx(T.domain);
            end
            if ( nargin < 2 )
                m = 1;
            end

            for l = 1:m % Refine m times.
                nDom = size(T.domain, 1);
                dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
                for k = 1:nDom
                    d = dom{k}; % Subdivide this square into two new pieces:
                    midPt  = mean(d(1:2));
                    dom{k} = [ d(1), midPt, d(3:4) ;
                               midPt, d(2), d(3:4) ];
                end
                T.domain = cell2mat(dom);              % Revert to a matrix.
                hMerge = reshape(1:2*nDom, 2, nDom).'; % New horizontal merge.
                T.mergeIdx = [hMerge, T.mergeIdx];     % Append to existing.
            end

        end

        function T = refiney(T, m)
        %REFINEY   Refine a domain in the y-direction.
        %   REFINEY(T) will divide each subdomain of T vertically into two
        %   new equally-sized pieces. The tree index information in the
        %   result is updated to reflect the new subdomains, which are the
        %   first to be merged.
        %
        %   REFINEY(T, M) will refine M times.

            if ( isempty(T.domain) )
                return
            end
            if ( isempty(T.mergeIdx) && length(T) > 1)
                warning('Empty tree index encountered. Building a default one.');
                T.mergeIdx = defaultIdx(T.domain);
            end
            if ( nargin < 2 )
                m = 1;
            end

            for l = 1:m % Refine m times.
                nDom = size(T.domain, 1);
                dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
                for k = 1:nDom
                    d = dom{k}; % Subdivide this square into two new pieces:
                    midPt  = mean(d(3:4));
                    dom{k} = [ d(1:2), d(3), midPt ;
                               d(1:2), midPt, d(4) ];
                end
                T.domain = cell2mat(dom);              % Revert to a matrix.
                hMerge = reshape(1:2*nDom, 2, nDom).'; % New horizontal merge.
                T.mergeIdx = [hMerge, T.mergeIdx];     % Append to existing.
            end

        end

        function T = removePatch(T, k)
        %REMOVEPATCH   Remove a patch (or patches) from a domain.
        %   T = REMOVEPATCH(T, K) removes the Kth patch(es) from the
        %   ultraSEMDomain T and updates the first level mergeIdx
        %   appropriately (by introducing NaNs). K may be a vector.

            % Sort for conveience:
            k = unique(sort(k(:)));
            % Remove specified patches:
            T.domain(k,:) = [];

            if ( size(T.domain, 1) == 0 )
                % We have removed all the patches!
                T = ultraSEMDomain();
                return
            end

            % Remove the corresponding indexes from the first level:
            idx1 = T.mergeIdx{1};
            for j = 1:numel(k)
                idx1(idx1 == k(j)) = NaN; % Replace missing patches with NaNs.
                idx1(idx1 > k(j)) = idx1(idx1 > k(j)) - 1; % Update.
                k = k - 1;
            end
            T.mergeIdx{1} = idx1;

        end

        function T = rot90(T, k)
        %ROT90   Rotate an ultraSEMDomain (clockwise) by 90 degrees.
        %   T = ROT90(T) rotates the ultraSEMDomain T clockwise by 90
        %   degrees.
        %
        %   T = ROT90(T, K) rotates by K*90 degrees.

            if ( nargin == 1 )
                k = 1;
            end

            if ( ~isinteger(k) )
                error('ULTRASEM:ULTRASEMDOMAIN:rot90:invalid', ...
                    'K must be an integer.');
            end

            k = mod(k, 4);
            if ( k == 1 )
                T.domain = T.domain(:, [3 4 1 2]);
                T = fliplr(T);
            elseif ( k == 2 )
                T = fliplr(flipud(T)); %#ok<FLUDLR>
            elseif ( k == 3 )
                T.domain = T.domain(:, [3 4 1 2]);
                T = flipud(T);
            else
                error('ULTRASEM:ULTRASEMDOMAIN:rot90:invalid', ...
                    'K must be an integer.');
            end

        end

        function T = transpose(T)
        %TRANSPOSE   Rotate an ultraSEMDomain clockwise by 90 degrees.
            T = rot90(T, -1);
        end

        function varargout = times(varargin)
        %TIMES   Scale an ultraSEMDomain.
        %   See MTIMES() for documentation.
        %
        % See also MTIMES.
            [varargout{1:nargout}] = mtimes(varargin{:});
        end

        function T = uplus( T )
        %UPLUS   Unary plus. Has no effect on an ultraSEMDomain.
        end

        function T = uminus( T )
        %UMINUS   Negate an ultraSEMDomain.
            T.domain = -T.domain;
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        % Construct alphabet domains:
        T = alphabet(str);

        function mergeIdx = defaultIdx( dom )
        %DEFAULTIDX   Make a default (naive) merge index for a domain.
        %   MERGEIDX = DEFAULTIDX(DOM) automatically constructs a merge
        %   index MERGEIDX for the domain DOM. The merge index will take
        %   the form
        %       MERGEIDX = {[1 2 ; 3 4 ; ... ; N-1 N], [1 2 ; ... ; N/2-1 N/2], ..., [1 2]},
        %   where N is the number of patches in DOM. NaNs are introduced
        %   where there are an odd number of patches on a level.
        %
        %   MERGEIDX = DEFAULTIDX(N) is equivalent.

            % Parse input:
            if ( isa(dom, 'ultraSEMDomain') )
                dom = dom.domain;
            end
            if ( isa(dom, 'scalar') )
                np = dom;             % Number of patches.
            else
                np = size(dom, 1);    % Number of patches.
            end
            nl = ceil(log2(np));      % Number of layers.
            mergeIdx = cell(1, nl);   % Initialise.

            for l = 1:nl
                if ( mod( np, 2 ) )
                    ii = [1:np, NaN]; % Pad with a NaN.
                else
                    ii = 1:np;
                end
                np = ceil( np/2 );    % Number of patches is halved each layer.
                mergeIdx{1,l} = reshape(ii, 2, np).';
            end

        end

        % Construct rectangular domains:
        T = rectangle(varargin);

        % Construct quadrilateral domains:
        T = quad(varargin);

        % Construct triangular domains:
        T = triangle(varargin);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
