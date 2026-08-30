function [R,Xsub,idx,E]=licols(X,tol)
% function [R,Xsub,idx,E]=licols(X,tol)
% Modified from licols, version 1.0.3, by Matt J.
% Original copyright (c) 2020, Matt J. Distributed under the BSD 3-Clause
% license; see LICENSE-licols.txt and THIRD_PARTY_NOTICES.md.
%Extract a linearly independent set of columns of a given matrix X
%INPUT:
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%OUTPUT:
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns; Xsub=X(:,idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Oct 30/21:  [~, R, E] = qr(full(X),0); % full(X) is now used
%  to guarantee diag(abs(R)) is nonincreasing
if ~nnz(X) %X has no non-zeros and hence no independent columns
    Xsub=[]; idx=[];R = [];
    return
end
if nargin<2, tol=1e-8; end   % 1e-8 used in order to get well-cond Xsub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, R, E] = qr(X,0);   % original way
% [~, R, E] = qr(full(X),0);   % full(X) now used to get nondecr diagr
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);
end
%Rank estimation  
rinds = diagr >= tol*max(diagr);
r = sum(rinds);
idx = sort(E(rinds));
Xsub=X(:,idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%  of function licols
