function v = mysvec(M)



n = length(M);
M(tril(true(n),-1)) = sqrt(2)*M(tril(true(n),-1));
v = M(tril(true(n)));

end