function [x,status_sub] = gsolve_fea(prob,u,d,modelsense,pars)
% feasibility version
% status_sub can take the following values
% 'SUCCESS'     we find a new affinely independent solution
% 'INFEASIBLE'  no affinely independent solution exists for this subproblem
% 'TIME_LIMIT'  MILP can't finish within time limit
% 'UNKNOWN'     solves failed for unknown reasons


tol = 1e-1;

% initialize output
% a new affinely indpendent vector if we found it
x = [];


% we try to find a feasible solution such that u^{T}x neq d
% we make this a feasibility problem by adding an extra constraint
% here, we preallocate the space for the extra constraint
prob.A = [prob.A; u'];
prob.rhs = [prob.rhs; d];
prob.obj = zeros(length(prob.obj),1);

% adjust the constraint
% we modify the constraints in two different ways
% i = 1, find u^{T}x >= d + tol
% i = 2, find u^{T}x <= d - tol
if strcmp(modelsense,'max')
    prob.rhs(end) = d+tol;
    prob.sense = [prob.sense; '>'];
else
    prob.rhs(end) = d-tol;
    prob.sense = [prob.sense; '<'];
end


% solver setting
gurobi_pars = [];
gurobi_pars.outputflag = 0;
gurobi_pars.TimeLimit = pars.sub_time;

% solve the feasibility problem
output = gurobi(prob,gurobi_pars);

% infeasible means it does not exist
if strcmp(output.status,'INFEASIBLE')
    status_sub = 'INFEASIBLE';
    return
end

% out of time
if strcmp(output.status,'TIME_LIMIT')
    status_sub = 'TIME_LIMIT';
    return
end


if ~strcmp(output.status,'OPTIMAL')
    warning('gsolve_fea:UnexpectedGurobiStatus', ...
        'Gurobi returned unexpected status: %s.', output.status);
    status_sub = 'UNKNOWN';
    return
end


% if the optimal value is indeed different from d
if abs(u'*output.x - d)>=tol
    % if abs(d) > 1e5
    % keyboard
    % end
    status_sub = 'SUCCESS';
    x = output.x;
    return
else
    % the solver should have reported infeasible
    status_sub = 'INFEASIBLE';
    x = output.x;
    return
end



end
