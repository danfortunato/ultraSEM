function [R, mergeIdx] = refine(R, m)
%REFINE   Refine ULTRASEM.RECT objects.    

mergeIdx = {};
% Parse inputs:
if ( nargin == 1 )
    m = 1;
elseif ( m == 0 )
    return
end

v = rectVertices(R);        % Get 1x4 vector of vertices.

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
    % Construct new merge indicies:
    hIdx = reshape(1:4*n, 2, 2*n).';   % New horizontal merge.
    vIdx = reshape(1:2*n, 2, n).';     % New vertical merge.
    mergeIdx = [hIdx, vIdx, mergeIdx]; %#ok<AGROW> Append to existing.  
end

R = ultraSEM.Rect(v);       % Assign new vertices..

end
