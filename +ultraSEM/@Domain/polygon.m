function T = polygon(v)
%ULTRASEM.POLYGON   Return a convex polygonal domain formed of quadrilaterals.
%   T = ULTRASEM.POLYGON(V) returns a convex polygonal ULTRASEM.DOMAIN T
%   with vertices V formed of N quads, where N is the polygon degree. The
%   vertices V should be given in anticlockwise order. If they are not,
%   they will be modified to be as such.
%
%   T = ULTRASEM.POLYGON(N) returns an ULTRASEM.DOMAIN T that is an N-sided
%   regular polygon.
%
%   It is important to note that each side of the polygon will be formed
%   from two adjacent grids.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    if ( isscalar(v) )
        % Construct a v-sided regular polygon.
        poly = nsidedpoly(v);
        T = ultraSEM.Domain.polygon(poly.Vertices);
        return
    end

    % Number of sides of the polygon
    n = size(v,1);

    % Check that the given points form a convex polygon:
    if ( ~all(unique(convhull(v)) == (1:n)') )
        error('ULTRASEM:DOMAIN:polygon:nonconvex', ...
            'Polygon is not convex.');
    end

    if ( ultraSEM.Domain.isClockwise(v) )
        % Switch the vertices to make it anticlockwise.
        v = flipud(v);
    end

    % Locate the center of the polygon:
    c = mean( v );

    % Construct quads:
    Q(n,1) = ultraSEM.Quad();
    for k = 1:n
        next = mod(k,n)+1;
        prev = mod(k-2,n)+1;
        vk = [ v(k,:) ; mean(v([k,next],:)) ; c ; mean(v([prev,k],:)) ];
        Q(k) = ultraSEM.Quad( vk );
    end

    % Construct merges:
    lvls = ceil(log2(n));
    mergeIdx = cell(1, lvls);
    for k = 1:lvls
        m = ceil(2^(1-k)*n);
        parray = padarray((1:m)', mod(m,2), NaN, 'post');
        mergeIdx{k} = reshape( parray, 2, ceil(m/2))';
    end

    % Construct domain:
    T = ultraSEM.Domain(Q, mergeIdx);

end

function b = padarray(varargin)
%PADARRAY Pad array.
%   B = PADARRAY(A,PADSIZE) pads array A with PADSIZE(k) number of zeros
%   along the k-th dimension of A.  PADSIZE should be a vector of
%   nonnegative integers.
%
%   B = PADARRAY(A,PADSIZE,PADVAL) pads array A with PADVAL (a scalar)
%   instead of with zeros.
%
%   B = PADARRAY(A,PADSIZE,PADVAL,DIRECTION) pads A in the direction
%   specified by DIRECTION. DIRECTION can be one of the following:
%
%       String or character vector values for DIRECTION
%       'pre'         Pads before the first array element along each
%                     dimension .
%       'post'        Pads after the last array element along each
%                     dimension.
%       'both'        Pads before the first array element and after the
%                     last array element along each dimension.
%
%   By default, DIRECTION is 'both'.
%
%   B = PADARRAY(A,PADSIZE,METHOD,DIRECTION) pads array A using the
%   specified METHOD.  METHOD can be one of the following:
%
%       String or character vector values for METHOD
%       'circular'    Pads with circular repetition of elements.
%       'replicate'   Repeats border elements of A.
%       'symmetric'   Pads array with mirror reflections of itself.
%
%   Class Support
%   -------------
%   When padding with a constant value, A can be numeric or logical.
%   When padding using the 'circular', 'replicate', or 'symmetric'
%   methods, A can be of any class.  B is of the same class as A.
%
%   Example
%   -------
%   Add three elements of padding to the beginning of a vector.  The
%   padding elements contain mirror copies of the array.
%
%       b = padarray([1 2 3 4],3,'symmetric','pre')
%
%   Add three elements of padding to the end of the first dimension of
%   the array and two elements of padding to the end of the second
%   dimension.  Use the value of the last array element as the padding
%   value.
%
%       B = padarray([1 2; 3 4],[3 2],'replicate','post')
%
%   Add three elements of padding to each dimension of a
%   three-dimensional array.  Each pad element contains the value 0.
%
%       A = [1 2; 3 4];
%       B = [5 6; 7 8];
%       C = cat(3,A,B)
%       D = padarray(C,[3 3],0,'both')
%
%   See also CIRCSHIFT, IMFILTER.

%   Copyright 1993-2017 The MathWorks, Inc.

args = matlab.images.internal.stringToChar(varargin);
[a, method, padSize, padVal, direction] = ParseInputs(args{:});

b = padarray_algo(a, padSize, method, padVal, direction);

end

%%%
%%% ParseInputs
%%%
function [a, method, padSize, padVal, direction] = ParseInputs(varargin)

narginchk(2,4);

% fixed syntax args
a         = varargin{1};
padSize   = varargin{2};

% default values
method    = 'constant';
padVal    = 0;
direction = 'both';

validateattributes(padSize, {'double'}, {'real' 'vector' 'nonnan' 'nonnegative' ...
    'integer'}, mfilename, 'PADSIZE', 2);

% Preprocess the padding size
if (numel(padSize) < ndims(a))
    padSize           = padSize(:);
    padSize(ndims(a)) = 0;
end

if nargin > 2
    
    firstStringToProcess = 3;
    
    if ~ischar(varargin{3})
        % Third input must be pad value.
        padVal = varargin{3};
        validateattributes(padVal, {'numeric' 'logical'}, {'scalar'}, ...
            mfilename, 'PADVAL', 3);
        
        firstStringToProcess = 4;
        
    end
    
    for k = firstStringToProcess:nargin
        validStrings = {'circular' 'replicate' 'symmetric' 'pre' ...
            'post' 'both'};
        string = validatestring(varargin{k}, validStrings, mfilename, ...
            'METHOD or DIRECTION', k);
        switch string
            case {'circular' 'replicate' 'symmetric'}
                method = string;
                
            case {'pre' 'post' 'both'}
                direction = string;
                
            otherwise
                error(message('images:padarray:unexpectedError'))
        end
    end
end

% Check the input array type
if strcmp(method,'constant') && ~(isnumeric(a) || islogical(a))
    error(message('images:padarray:badTypeForConstantPadding'))
end

end

function b = padarray_algo(a, padSize, method, padVal, direction)
%PADARRAY_ALGO Pad array.
%   B = PADARRAY_AGLO(A,PADSIZE,METHOD,PADVAL,DIRECTION) internal helper
%   function for PADARRAY, which performs no input validation.  See the
%   help for PADARRAY for the description of input arguments, class
%   support, and examples.

%   Copyright 2014 The MathWorks, Inc.

if isempty(a)
    
    numDims = numel(padSize);
    sizeB = zeros(1,numDims);
    
    for k = 1: numDims
        % treat empty matrix similar for any method
        if strcmp(direction,'both')
            sizeB(k) = size(a,k) + 2*padSize(k);
        else
            sizeB(k) = size(a,k) + padSize(k);
        end
    end
    
    b = mkconstarray(class(a), padVal, sizeB);
    
elseif strcmpi(method,'constant')
    
    % constant value padding with padVal
    b = ConstantPad(a, padSize, padVal, direction);
else
    
    % compute indices then index into input image
    aSize = size(a);
    aIdx = getPaddingIndices(aSize,padSize,method,direction);
    b = a(aIdx{:});
end

if islogical(a)
    b = logical(b);
end

end

%%%
%%% ConstantPad
%%%
function b = ConstantPad(a, padSize, padVal, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
            
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
            
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;

end

function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2013 The MathWorks, Inc.  

out = repmat(feval(class, value), size);

end
