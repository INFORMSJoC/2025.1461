function output_warmup = find_initial_fea(prob,pars,Vbar)
% obtain an initial feasible solution for prob

% get the maximum number of affinely indepedent solutions
if strcmp(pars.var_mode,'B')
    K = sum(pars.bvar) + 1 - size(Vbar,2);
else
    K = length(prob.obj) + 1 - size(Vbar,2);
end

% the number of variables
n = length(prob.obj);

% change the problem into a feasibility problem
prob.obj = zeros(n,1);

% get single feasible solution
[V,status] = get_one(prob,pars);
nr_sub = 1;

% if warm_up is true
% then we obtain some more zero-one solutions
if and(pars.warm_up == true,~isempty(V))
    % if V(i,:) is all-zeros or all-ones,
    for i = 1:sum(get_bvar(prob))
        for j = 0:1
            % if xi is all j
            if all(V(i,:) == j)
                % try to find a solution so that xi = 1-j
                myprint(pars.DEBUG,'finding sol so that %d-th entry is %d\n',i,1-j)
                V = get_xi(prob,i,1-j,V,pars);
                nr_sub = nr_sub + 1;
            end
        end
    end
end

% save the results
r0 = size(V,2); % number of feasible solution from warmup
output_warmup.r0 = r0;
output_warmup.V = V;
output_warmup.nr_sub = nr_sub; % number of subroblems solved from warmup
output_warmup.size = size(V,2);

% obtain the status
if strcmp(status,'INFEASIBLE')
    % infeasible problem
    output_warmup.status = 'INFEA';
end

if r0 == K
    % already found affine
    output_warmup.status = 'SUCC_SIZE';
end

if strcmp(status,'TIME_LIMIT')
    % reached time limit
    output_warmup.status = status;
end

myprint(pars.DEBUG,'warmup yields %d feasible solutions\n',size(V,2))
% output.time = toc(total_time);

end

function V = get_xi(prob,i,j,V,pars)

% the number of variables
n = length(prob.obj);

prob.rhs = prob.rhs - prob.A(:,i)*j;
prob.A(:,i) = [];
prob.vtype(i) = [];
prob.lb(i) = [];
prob.ub(i) = [];
prob.obj(i) = [];

% solve the problem
g_pars = [];
g_pars.outputflag = 0;
g_pars.TimeLimit = pars.sub_time;
g_output = gurobi(prob,g_pars);

status = g_output.status;

if strcmp(status,'OPTIMAL')
    x = zeros(n,1);
    x(1:i-1) = g_output.x(1:i-1);
    x(i+1:end) = g_output.x(i:end);
    x(i) = j;
    V = [V x];
    % elseif strcmp(status,'INFEASIBLE')
    % fprintf('there does not exist feasible x such that x(%d) = %d\n',i,j)
end


end


function [V,status] = get_one(prob,pars)

% the number of variables
n = length(prob.obj);

% change the problem into a feasibility problem
prob.obj = zeros(n,1);

% solve the problem
g_pars = [];
g_pars.outputflag = 0;
g_pars.TimeLimit = pars.sub_time;
g_out = gurobi(prob,g_pars);
status = g_out.status;
% save the result
if strcmp(status,'OPTIMAL')
    myprint(pars.DEBUG,'An initial solution is obtained\n')
    V = g_out.x; % save the feasible solution
else
    myprint(pars.DEBUG,'Failed to find an initial solution\n')
    V = [];
end

end