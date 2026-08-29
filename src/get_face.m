function [Vnew,d,u,r,U] = get_face(V,pars,tol)
% Vnew = the range space of V
% d = the smallest positive singular
% u = the left singular vector associated with d
% r = number of positive eigenvalues

Vnew = [];
d = [];
u = [];
r = [];

if isempty(V)
    return
end


% only for binary variables or all variables
if strcmp(pars.var_mode,'B')
    bvar = pars.bvar;
else
    bvar = 1:size(V,1);
end

% get the face
X = V(bvar,:);
k = size(X,2);
X = [ones(1,k);X]/sqrt(k);

% get the svd
[P,S,~] = svd(X,'econ');
if nargin == 3
    tol = max(tol,max(size(X)) * eps(S(1,1)));
else
    tol = max(size(X)) * eps(S(1,1));
end
Vnew = P(:, diag(S) > tol);

% get the minimum singular value
[d,idx] = min(diag(S));
u = P(:,idx);
U = P(:,diag(S)<1e-1);
r = size(Vnew,2);

% if full rank, then choose Vb to be the identity
if size(Vnew,1) == size(Vnew,2)
    Vnew = eye(size(Vnew,1));
end


end
