function output = get_feasible_sols(prob,pars)
% given a MILP problem in gurobi format, this function generates a maximal
% affinely independent set of of feasible solutions saved in V
% output.status
% 'SUCCESS'    - V spans the affine hull
% 'TIME_LIMIT' - MILP does not solve the problem within TIME_LIMIT
% 'NUMERIC'    - MILP does not return correct solution somehow
% 'UNKNOWN'    - MILP failed to find a solution for unknown reasons

myprint(pars.DEBUG,'\nMILP FR started\n')

% load the settings
pars = load_setting(prob,pars);

% get the implicit equalities
[Vbar,output_implict] = get_implicit_equalities(prob,pars);

% restrict the feasible set to a well-conditioned range.
% prefer numerical stability over reduction size
% note that this only changes the FR, not the original problem.
n = length(prob.lb);
prob.ub = min(prob.ub,1e3*ones(n,1));
prob.lb = max(prob.lb,-1e3*ones(n,1));

% obtain an initial set of affinely independent feasible solutions
output_warmup = find_initial_fea(prob,pars,Vbar);

% find a set of affinely independent feasible solutions
output_aff = find_V(prob,pars,output_warmup,Vbar);

% add more feasible solutions so that the smallest eigenvalue is increased
output_refine = refine(prob,pars,output_aff);

% save the output
output.implicit = output_implict;
output.warmup = output_warmup;
output.aff = output_aff;
output.refine = output_refine;
% keyboard
end

function output = find_V(prob,pars,output_warmup,Vbar)

total_time = tic;

% get the maximum number of affinely indepedent solutions
if strcmp(pars.var_mode,'B')
    K = sum(pars.bvar) + 1 - size(Vbar,2);
else
    K = length(prob.obj) + 1 - size(Vbar,2);
end

V = output_warmup.V;
output.nr_sub = output_warmup.nr_sub;

% if output has field `status', then we terminate as one of the following
% three states [SUCCESS,TIME_LIMIT,INFEASIBLE]
% otherwise we proceed
if ~isfield(output_warmup,'status')
    % at each iteration, we try to find a new affinely independent solution
    % at the end of each iteration, we have three possible states
    % 1. we found a new solution
    % 2. we prove that the current V spans the affine hull
    % 3. the algorithm failed to generate a new solution
    myprint(pars.DEBUG,'\n',[])
    r0 = size(V,2);

    % if warmup already yields a maximal number of feasible solutions
    if r0+1 > K
        myprint(pars.DEBUG,'Maximality proved by size.\n')
        output.status = 'SUCC_SIZE';
    else
        % generate more feasible solutions
        for i = r0+1:K
            % print the title
            myprint(pars.DEBUG,'Iteration i = %d\n',i-r0)
            myprint(pars.DEBUG,'    current size = %d, max size = %d\n',size(V,2),K)
            myprint(pars.DEBUG,'    current time = %.2f, max time = %d\n',toc(total_time),pars.max_time)


            % try to find a new affinely independent vector
            [V,status_iter,nr_sub] = find_new_sol(prob,V,Vbar,pars,total_time);

            % remember how many subproblems are solved
            output.nr_sub = output.nr_sub + nr_sub;
            output.V = V;

            % check the stopping condition
            if and(strcmp(status_iter,'SUCCESS'),size(V,2) == K)
                % we find n + 1 affinely independent sol.
                myprint(pars.DEBUG,'Maximality proved by size.\n')
                output.status = 'SUCC_SIZE';
            elseif strcmp(status_iter,'MAXIMAL')
                % MILP proves that V is a maximally affinely independent solutions
                myprint(pars.DEBUG,'Maximality proved by MILP\n')
                output.status = 'SUCC_ILP';
            elseif any(strcmp(status_iter,{'TIME_LIMIT', 'NUMERIC', 'UNKNOWN'}))
                % the algorithm fails
                myprint(pars.DEBUG,'Failed to generate due to %s\n',status_iter)
                output.status = status_iter;
            elseif toc(total_time) > pars.max_time
                myprint(pars.DEBUG,'Reached global time limit\n')
                output.status = 'TIME_LIMIT';
            end

            % if output status has a value, then we terminate
            if isfield(output,'status')
                break
            end
        end
    end
else
    % if warmup already has a status, then set output.status the same
    output.status = output_warmup.status;
end

% save the results
output.V = V;
output.size = size(V,2); % number of feasible solutions
[~,d,~,rk] = get_face(V,pars);
output.rk = rk; % rank of the feasible solutions
output.cond = d; % the smallest positive singular value
output.time = toc(total_time); % total time

end



function [V,status_iter,nr_sub] = find_new_sol(prob,V,Vbar,pars,total_time)
% output:
%        V           := a new set of affinely independent feasible solutions
%        status_iter := the status of the current iteration
%        nr_sub      := number of MILP problems being solved at this iteration
%

% generate a new solution
[V,status_iter,nr_sub] = find_x(prob,V,Vbar,pars);

% if we couldn't find a new solution within time limit, and we are using
% the random objective and within global time limit
% then we give one more try
if strcmp(status_iter,'TIME_LIMIT')&&pars.rand&&(toc(total_time) < pars.max_time)
    [V,status_iter,nr_sub1] = find_x(prob,V,Vbar,pars);
    nr_sub = nr_sub + nr_sub1;
    myprint(pars.DEBUG,'       generate a different random objective\n')
    myprint(pars.DEBUG,'       current time = %.2f, max time = %d\n',toc(total_time),pars.max_time)
end

end

function [V,status_iter,nr_sub] = find_x(prob,V,Vbar,pars)
% status_iter is the status of the current iteration for generating a new
% affinely independent solution, and it can be the following values:
% 'SUCCESS'    = find a new affinely independent solution
% 'MAXIMAL'    = V is proved to be maximal affinely independent solutions
% 'TIME_LIMIT' = MILP failed to generate a new solution within time limit
% 'NUMERIC'    = MILP generates a new solution, but somehow it is incorrect
% 'UNKNOWN'    = MILP failed to find a solution for unknown reasons


% parameters
DEBUG = pars.DEBUG;


% generate the search directions
if strcmp(pars.var_mode,'B')
    U = find_orth(V,Vbar,pars.bvar);     % compute the linear subspace U
else
    U = find_orth(V,Vbar);               % compute the linear subspace U
end
if pars.rand == true
    % we only consider a random linearly combination
    k = size(U,2);          % the number of columns in U
    U = U*randi([1 k],k,1); % generate a random linear combination of U
end
k = size(U,2);      % the number of columns in U
d = U'*V;           % compute the constants d
nr_sub = 0;         % total number of subproblems being solved



myprint(DEBUG,'       the orthogonality is %.4e\n',norm(U'*(V(:,1:end-1) - V(:,end)),'fro') )
myprint(DEBUG,'       #subproblems = %d max sub_time = %d\n',k,pars.sub_time)


% try to find a new affinely independent sol within the time limit
% for each column in U
for i = 1:k
    % for max/min problems
    for j = 1:2
        % update total number of subproblems being solved
        nr_sub = nr_sub + 1;

        % set the model sense
        if j == 1
            modelsense = 'max';
        else
            modelsense = 'min';
        end

        myprint(DEBUG,'       the %d-th %s problem with target %8.5f ',i,modelsense,d(i))

        % solve the problem based on the mode
        switch pars.mode
            case 'opt'
                [x,status_sub] = gsolve_opt(prob,U(:,i),d(i),modelsense,pars);
            case 'fea'
                [x,status_sub] = gsolve_fea(prob,U(:,i),d(i),modelsense,pars);
        end


        % decide the results
        if strcmp(status_sub,'SUCCESS')
            % if we found a solution
            myprint(DEBUG,'(OPTIMAL) - (obj %.4f) - (diff %.4f)\n',U(:,i)'*x,abs(U(:,i)'*x-d(i)))
            % if norm(x) > 1e5
            %     keyboard
            % end
            V = [V x];
            status_iter = 'SUCCESS';
            return
        elseif strcmp(status_sub,'INFEASIBLE')
            % the problem does not contain a solution so that u'(x-v0) < or > 0
            myprint(DEBUG,'(INFEASIBLE)\n')
        elseif strcmp(status_sub,'TIME_LIMIT')
            % if time limit is reached,
            status_iter = 'TIME_LIMIT';
            myprint(DEBUG,'(TIME_LIMT)\n')
            return
        elseif strcmp(status_sub,'SCALE')
            myprint(DEBUG,'(SCALE)\n')
            status_iter = 'SCALE';
            return
        elseif strcmp(status_sub,'UNKNOWN')
            myprint(DEBUG,'(UNKNOWN)\n')
            status_iter = 'UNKNOWN';
            return
        else
            error('get_feasible_sols:UnexpectedSubproblemStatus', ...
                'Unexpected subproblem status: %s.', status_sub);
        end
    end
end


% at this point, we know all subproblems are INFEASIBLE
myprint(DEBUG,'       ALL subproblems are INFEASIBLE\n')
% then we prove that V is maximal affinely independent
status_iter = 'MAXIMAL';

end


function output = refine(prob,pars,output_aff)
% prob := the problem
% V    := a given set of feasible solutions
% t1   := the time has been used already so far
% pars := parameters

total_time = tic;
t1 = output_aff.time;
V = output_aff.V;
% if V has 1 or less columns, then there is nothing to refine
if or(size(V,2)==0,pars.refine==false)
    output.V = -1;
    output.nr_sub = -1;
    output.size = -1;
    output.rk = -1;
    output.cond = -1;
    output.time = toc(total_time);
    return
elseif or(size(V,2) == 1,output_aff.cond > 1e-3)
    output = output_aff;
    output.nr_sub = 0;
    output.time = 0;
end

% switch between biary or all variable modes
if strcmp(pars.var_mode,'B')
    N = sum(or(prob.vtype=='I',prob.vtype=='B')); % number of binary varibles;
elseif strcmp(pars.var_mode,'A')
    N = length(prob.obj);
end

% solver setting
gurobi_pars = [];
gurobi_pars.outputflag = 0;
gurobi_pars.TimeLimit = 10;
gurobi_pars.PoolSolutions = N;

% maximization problem
prob.modelsense = 'max';

% get the smallest positive singularvalue
[~,d,u,rk,U1] = get_face(V,pars);


n = length(prob.obj);
nb = sum(or(prob.vtype=='I',prob.vtype=='B'));

% add more solutions
for nr_sub = 1:3*N
    % set up the objective
    if strcmp(pars.var_mode,'B')
        prob.obj = zeros(n,1);
        prob.obj(1:nb) = u(2:end);
    else
        prob.obj = u(2:end);
    end
    % prob.obj = sum(U1(2:end,:),2);

    % solve the optimization problem
    g_output = gurobi(prob,gurobi_pars);

    % get the solution
    xn = get_xn(g_output);
    if isempty(xn)
        break
    end

    % add the new solution
    V1 = [V xn(:,1)];
    % [d1,u1,rk1] = get_min_eig(V1,pars);
    [~,d1,u1,rk1,U1] = get_face(V1,pars,1e-8);
    % fprintf('\nrank %d new rank %d eig %.6e new eig %.6e\n',rk,rk1,d,d1)

    % if the rank is increased, or the smallest eigenvalue is improved
    % then we continue; otherwise, we terminate
    if or(and(rk==rk1,d1<=d),rk1<rk)
        break
    end

    % update the solution
    V = V1;
    d = d1;
    u = u1;
    rk = rk1;

    % the smallest singular value is at least 1e-3
    if d > 1e-3
        break
    end

    % the total running time is limited
    if toc(total_time) > pars.max_time - t1
        break
    end
end

% % extract linearly independent columns from V
% [~,S,~] = svd(V,'econ');
% tol = max(size(V)) * eps(S(1,1));
% [~,V] = licols(V,tol);
% [~,d,~,rk] = get_face(V,pars);

% save the output
output.V = V;
output.nr_sub = nr_sub;
output.size = size(V,2);
output.rk = rk;
output.cond = d;
output.time = toc(total_time);


end


function X = get_xn(g_out)
if ~isfield(g_out,'pool')
    X = [];
    return
end

pool = g_out.pool;

n = length(pool);                % Number of struct entries
vlen = length(pool(1).xn);       % Length of each vector
X = zeros(vlen, n);           % Preallocate matrix

for i = 1:n
    X(:, i) = pool(i).xn;        % Fill column i with S(i).xn
end

end



function pars = load_setting(prob,pars)
% get the default parameters
if ~isfield(pars,'DEBUG')
    pars.DEBUG = false;
end

if ~isfield(pars,'mode')
    pars.mode = 'opt';
end

if ~isfield(pars,'var_mode')
    pars.var_mode = 'B';
end

if ~isfield(pars,'rand')
    pars.rand = false;
end

if ~isfield(pars,'warm_up')
    pars.warm_up = false;
end

if ~isfield(pars,'solpool')
    pars.solpool = false;
end

if ~isfield(pars,'max_time')
    pars.max_time = false;
end

if ~isfield(pars,'sub_time')
    pars.sub_time = 5;
end

if ~isfield(pars,'trial')
    pars.trial = 1;
end

if strcmp(pars.var_mode,'B')
    pars.bvar = get_bvar(prob);
end

end
