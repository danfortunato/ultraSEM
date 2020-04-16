function chebpolyfill(sol)
%CHEBPOLYFILL   Display coefficients graphically.
%   CHEBPOLYFILL(SOL) plots the absolute values of the coefficients
%   underlying the representation of the ULTRASEM.SOL SOL on a semilog
%   scale using filled polygons.

holdState = ishold();
n = size(sol.u{1}, 1);

% Loop over the patches:
data = {};
for k = 1:length(sol)
    u = abs(sol.u{k});
    tail = [u(end-floor(.1*n):end,:) ; u(:,end-floor(.1*n):end)'];
    err(k) = norm(tail(:),inf)./max(u(:));
    err(k) = log10(err(k));
    d = sol.domain(k,:);
    data = [data, d([1 1 2 2 1]), d([3 4 4 3 3]), err(k)];
end

fill(data{:});

if ( ~holdState )
    hold off;
end

end
