function U = find_orth(V,Vbar,bvar)
% find the orthogonal complement of V(:,1:r-1) - V(:,r)
% the columns of U span the orthogonal complement of V(1:r-1)-V(:,r)


if nargin == 3
    % if we only consider binary variables that are specified by Bidx
    Vnew = [V(bvar,:) Vbar];
    Unew = my_orth(Vnew);
    U = zeros(size(V,1),size(Unew,2));
    U(bvar,:) = Unew;
else
    % if we consider all variables
    U = my_orth([V Vbar]);
end

% remarks:
%     U = null(L','r');
% with 'r'; it uses rank-revealing Q, leading to very inaccurate solution
% Althought [V U] still spans Rn, but we really need V'*U = 0 accurately
% thus, this command should be avoided.

end


function U = my_orth(V)
[n,r] = size(V);
if r == 1
    % identity
    U = eye(n);
else
    % compute the linear subspace L
    L = V(:,1:r-1) - V(:,r);
    % compute U
    U = null(full(L)');
end
return

end