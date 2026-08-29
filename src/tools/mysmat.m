function M = mysmat(v,n)

M = zeros(n);
M(tril(true(n))) = v;
M(tril(true(n),-1)) = M(tril(true(n),-1))/sqrt(2);
M = M + tril(M,-1)';

end