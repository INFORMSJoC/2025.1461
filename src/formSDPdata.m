function [A1,B1,b1,A2,B2,b2,nb] = formSDPdata(prob,pars)
% form an SDP relaxation for the mixed-binary non-convex quadratic program:
% the constraints are of the form
% <A,X> + By <=> b
% where sense is a vector of chars indicating the constraint is <, = or >

% convert the problem into Aineq*x<=bineq and Aeq*x = beq format
[Aineq,bineq,Aeq,beq] = FRAformat(prob);

% the number of binary and all variables
n = length(prob.obj);
nb = sum(or(prob.vtype=='B',prob.vtype=='I'));
nc = n - nb;

% number of equality and inequalities
m1 = size(Aeq,1);
m2 = size(Aineq,1);



switch pars.sdptype
    case 'SHOR'
        A1 = cell(2,1);
        b1 = cell(2,1);
        % get the equality constraints Aeq*x = beq
        [A1{1},b1{1}] = get_linear_constraints(Aeq,beq);

        % get the arrow constraints
        [A1{2},b1{2}] = get_arrow(prob);

        % get the inequality constraints
        [A2,b2] = get_linear_constraints(Aineq,bineq);

        A1 = cell2mat(A1);
        b1 = cell2mat(b1);
        B1 = [];
        B2 = [];


    case 'SHOR_B'
        A1 = cell(2,1);
        B1 = cell(2,1);
        b1 = cell(2,1);
        % get the linear equality constraints
        [A1{1},B1{1},b1{1}] = get_linear_constraints_B(Aeq,beq,nb);

        % get the arrow constraints
        [A1{2},B1{2},b1{2}] = get_arrow_B(n,nb);

        % get the linear inequality constraints
        [A2,B2,b2] = get_linear_constraints_B(Aineq,bineq,nb);

        A1 = cell2mat(A1);
        B1 = cell2mat(B1);
        b1 = cell2mat(b1);


    case 'DNN'
        A1 = cell(3,1);
        b1 = cell(3,1);
        A2 = cell(3,1);
        b2 = cell(3,1);
        % get the equality constraints Aeq*x = beq
        % tic
        [A1{1},b1{1}] = get_linear_constraints(Aeq,beq);
        % toc

        % get the equality constraints <aa',X> = b^2
        % tic
        [A1{2},b1{2}] = get_equality2(Aeq,beq);
        % toc

        % get the arrow constraints
        % tic
        [A1{3},b1{3}] = get_arrow(prob);
        % toc
    
        % get the inequality constraints
        % tic
        [A2{1},b2{1}] = get_linear_constraints(Aineq,bineq);
        % toc

        % [A2{2},b2{2}] = get_non_negative_B(prob);

        % get the RLT lower bounds
        % tic
        % [C2,d2] = get_lower_bounds(prob);
        % toc

        % get the RLT upper bounds
        % tic
        % [C3,d3] = get_upper_bounds(prob);
        % toc
        % keyboard

        % tic
        [A2{2},b2{2}] = get_rlt_bounds(prob,prob.lb);
        % toc

        % tic
        [A2{3},b2{3}] = get_rlt_bounds(prob,prob.ub);
        % toc
        % 
        % norm(A2{3} - C3,'fro')
        % norm(b2{3}-d3,'fro')
        % norm(A2{2} - C2,'fro')
        % norm(b2{2}-d2,'fro')
      
        % 
        A1 = cell2mat(A1);
        b1 = cell2mat(b1);
        B1 = [];
        A2 = cell2mat(A2);
        b2 = cell2mat(b2);
        B2 = [];
end
% keyboard
end



function [A1,B1,b1] = get_linear_constraints_B(A,b,nb)

m = size(A,1);
% the number of scalar variables in the matrix variable
n1 = nchoosek(nb+2,2);

% obtain the constraints in vector form
% the equality constraints
A1 = zeros(m,n1);
for i = 1:m
    a = A(i,1:nb);
    A1(i,:) = mysvec(sparse([ones(1,nb) 2:nb+1],[2:nb+1 ones(1,nb)],[a/2 a/2],nb+1,nb+1));
end
B1 = A(:,nb+1:end);
b1 = b;

end


function [A,B,b] = get_arrow_B(n,nb)

n1 = nchoosek(nb+2,2);

% the arrow constraints
A = zeros(nb+1,n1);
A(1,:) = mysvec(sparse(1,1,1,nb+1,nb+1));
for i = 1:nb
    A(i+1,:) = mysvec(sparse([i 1 i],[i i 1],[1 -1/2 -1/2],nb+1,nb+1));
end
B = zeros(nb+1,n-nb);
b = [1;zeros(nb,1)];

end





function [A1,b1] = get_linear_constraints(A,b)

[m,n] = size(A);
% the number of scalar variables in the matrix variable
n1 = nchoosek(n+2,2);

% obtain the constraints in vector form
% the equality constraints
A1 = zeros(m,n1);
for i = 1:m
    A1(i,:) = mysvec(sparse([ones(1,n) 2:n+1],[2:n+1 ones(1,n)],[A(i,:)/2 A(i,:)/2],n+1,n+1));
end
b1 = b;

end


function [A,b] = get_equality2(Aeq,beq)

[m,n] = size(Aeq);
% the number of scalar variables in the matrix variable
n1 = nchoosek(n+2,2);

% obtain the constraints in vector form
% the equality constraints
A = zeros(m,n1);
for i = 1:m
    a = Aeq(i,:);
    A(i,:) = mysvec([0;a']*[0 a]);
end
b = beq.^2;

end

function [A,b] = get_arrow(prob)

n = length(prob.obj);
nb = sum(or(prob.vtype == 'B',prob.vtype == 'I'));
n1 = nchoosek(n+2,2);

% the arrow constraints for the binary variables
A = zeros(nb+1,n1);
A(1,:) = mysvec(sparse(1,1,1,n+1,n+1));
for i = 2:nb+1
    A(i,:) = mysvec(sparse([i 1 i],[i i 1],[1 -1/2 -1/2],n+1,n+1));
end
b = [1;zeros(nb,1)];

end


function [A,b] = get_non_negative_B(prob)
% impose the non-negativity constraints

n = length(prob.obj);
nb = sum(or(prob.vtype == 'B',prob.vtype == 'I'));
n1 = nchoosek(n+2,2);
m1 = nchoosek(nb,2);

A = zeros(m1,n1);
b = zeros(m1,1);
k = 0;
for i = 2:nb+1
    for j = i+1:nb+1
            k = k + 1;
            A(k,:) = mysvec(sparse([i j],[i j],...
                                    [-1 -1],n+1,n+1));
            b(k) = 0;
    end
end

if length(b)~=k
    error('formSDPdata:NonnegativeConstraintCountMismatch', ...
        'Expected %d nonnegative constraints but constructed %d.', ...
        length(b), k);
end

end

function [A,b] = get_rlt_bounds(prob,bnd)
% impose the RLT type lower bounds
% xi - li >=0
% xj - lj > =0 
% yields
% xixj - xilj - xjli + lilj>=0
% or equivalenlty
% -xixj + xilj + xjli <= lilj
%
% if i = j, it becomes
% (-1 + li + li)xi <= li^2
%
% for i neq j, and li*lj finite
% each constraint has 
% (i,j)th entry -1/2,
% (j,i)th entry -1/2
% (0,i)th entry lj/2
% (i,0)th entry lj/2
% (0,j)th entry li/2
% (j,0)th entry li/2
% right hand side li*lj

% here n is the size of the matrix
n = length(prob.obj) + 1;
n1 = nchoosek(n+1,2); % the number of lower triangular elements 

% get the bounds
bnd = [Inf;bnd];

I = kron(ones(n,1),[1:n]');
J = kron([1:n]',ones(n,1));
bnd2 = [kron(bnd,ones(n,1)) kron(ones(n,1),bnd)];

% only multiply constraints i<j
idx = I >= J;
I = I(idx,:);
J = J(idx,:);
bnd2 = bnd2(idx,:);

% only keep constraints  li*lj finite
idx = and(isfinite(bnd2(:,1)),isfinite(bnd2(:,2)));
I = I(idx,:);
J = J(idx,:);
bnd2 = bnd2(idx,:);

m = length(I);

% get the product of bounds
b = bnd2(:,1).*bnd2(:,2);

% the coefficient x_{ij} is in IND_ij entry and has value IND_k 
IND_ij1 = lower_tri_index(I, J, n);
IND_k1 = ones(m,1);
IND_k1(I == J) = -1; % diagonal entries is -1
IND_k1(I ~= J) = sqrt(2)*(-1/2); % off-diagonal entries is sqrt(2)*(-1/2)

% the coefficient x_{i,1} is in IND_ij2 entry and has value IND_k2
IND_ij2 = lower_tri_index(I, ones(m,1), n);
IND_k2 = sqrt(2)*(1/2)*bnd(J);

% the coefficient x_{j,1} is in IND_ij2 entry and has value IND_k2
IND_ij3 = lower_tri_index(J, ones(m,1), n);
IND_k3 = sqrt(2)*(1/2)*bnd(I);

A = sparse(kron(ones(3,1),[1:m]'),[IND_ij1;IND_ij2;IND_ij3],[IND_k1;IND_k2;IND_k3],m,n1);


end


function [A,b] = get_lower_bounds(prob)
% impose the RLT type lower bounds
% xi - li >=0
% xj - lj > =0 
% yields
% xixj - xilj - xjli + lilj>=0
% or equivalenlty
% -xixj + xilj + xjli <= lilj


n = length(prob.obj);
n1 = nchoosek(n+2,2);

% get the lower bounds
lb = [Inf;prob.lb];
LB = lb*lb'; % get li*lj
m1 = (sum(sum(isfinite(LB))) - sum(isfinite(diag(LB))))/2 + sum(isfinite(diag(LB)));

A = zeros(m1,n1);
b = zeros(m1,1);
k = 0;
for i = 2:n+1
    for j = i:n+1
        if isfinite(LB(i,j))
            k = k + 1;
            A(k,:) = mysvec(sparse([i j 1 i 1 j],[j i i 1 j 1],...
                                    [-1/2 -1/2 lb(j)/2 lb(j)/2 lb(i)/2 lb(i)/2],n+1,n+1));
            b(k) = LB(i,j);
        end
    end
end

if length(b)~=k
    error('formSDPdata:LowerBoundConstraintCountMismatch', ...
        'Expected %d lower-bound RLT constraints but constructed %d.', ...
        length(b), k);
end

end



function [A,b] = get_upper_bounds(prob)
% impose the RLT type upper bounds
% ui - xi >= 0
% uj - xj >= 0 
% yields
% xi*xj - xi*uj - xj*ui + ui*uj>=0
% or equivalently
% -xi*xj + xi*uj + xj*ui <= ui*uj


n = length(prob.obj);
n1 = nchoosek(n+2,2);

% get the lower bounds
ub = [Inf;prob.ub];
UB = ub*ub'; % get li*lj
m1 = (sum(sum(isfinite(UB))) - sum(isfinite(diag(UB))))/2 + sum(isfinite(diag(UB)));

A = zeros(m1,n1);
b = zeros(m1,1);
k = 0;
for i = 2:n+1
    for j = i:n+1
        if isfinite(UB(i,j))
            k = k + 1;
            A(k,:) = mysvec(sparse([i j 1 i 1 j],[j i i 1 j 1],...
                                    [-1/2 -1/2 ub(j)/2 ub(j)/2 ub(i)/2 ub(i)/2],n+1,n+1));
            b(k) = UB(i,j);
        end
    end
end

if length(b)~=k
    error('formSDPdata:UpperBoundConstraintCountMismatch', ...
        'Expected %d upper-bound RLT constraints but constructed %d.', ...
        length(b), k);
end


end



% function [Anew,Bnew,bneq] = obtain_RLT(Anew,b,nb)
%
% [m,n] = size(Anew);
% nc = n - nb;
%
% Anew = cell(m*(nb+1),1);
% Bnew = zeros(m*(nb+1),nc*(nb+1));
% bnew = zeros(m*(nb+1),1);
% k = 0;
% for i = 1:m
%     ab = Anew(i,1:nb);
%     ac = Anew(i,nb+1:end);
%     for j = 0:nb
%         k = k + 1;
%         temp = zeros(nb+1);
%         if j == 0
%             temp(j+1,:) = [b(i) ab];
%         else
%             temp(j+1,:) = [0 ab];
%             temp(j+1,j+1) = temp(j+1,j+1) - b(j);
%         end
%         temp = (temp +temp')/2;
%         Anew{k} = temp;
%         Bnew(k,:) = kron(sparse(1,j+1,1,1,nb+1),ac);
%     end
% end
%
%
% end

% % obtain the equality constraints in matrix form
% A1 = cell(m1 + nb + 1,1);
% for i = 1:m1
%     A1{i} = sparse(2:nb+1,2:nb+1,Aeq(i,1:nb),nb+1,nb+1);
% end
% B1 = [Aeq(:,nb+1:end);zeros(nb+1,n-nb)];
% b1 = [beq;1;zeros(nb,1)];
% A1{m1+1} = sparse(1,1,1,nb+1,nb+1);
% for i = 1:nb
%     A1{m1+1+i} = sparse([i 1 i],[i i 1],[1 -1/2 -1/2],nb+1,nb+1);
% end
%
% A2 = cell(m2,1);
% for i = 1:m2
%     A2{i} = sparse(2:nb+1,2:nb+1,Aineq(i,1:nb),nb+1,nb+1);
% end
% B2 = Aineq(:,nb+1:end);
% b2 = bineq;
