function output = solve_SDP(prob,pars)

% get the SDP relaxation
% [A1,B1,b1,A2,B2,b2,nb] = get_SDPrelaxation(prob,pars);

[A1,B1,b1,A2,B2,b2] = formSDPdata(prob,pars);

n = length(prob.obj);

% apply FR
if pars.FR
    % % obtain the range space
    % W = pars.W;
    % V = null(W);
    % r = size(V,2);
    %
    % A1 = compute_VAV(A1,V);
    % A2 = compute_VAV(A2,V);
    %
    % X = sdpvar(r);
    %
    % [~,Asub] = licols([A1 b1]');
    % Asub = Asub';
    % A1 = Asub(:,1:end-1);
    % b1 = Asub(:,end);
    %
    %
    % % remove tiny entries after FR reformulation
    % A1(abs(A1)<1e-12) = 0;
    % A2(abs(A2)<1e-12) = 0;
    % b1(abs(b1)<1e-12) = 0;
    % b2(abs(b2)<1e-12) = 0;
    %
    %
    %
    % % set up the model
    % F = [X>=0; A1*mysvec(X) == b1; A2*mysvec(X) <= b2];
    %
    % % the objective function
    % obj = mysvec(V'*blkdiag(0,prob.Q)*V)'*mysvec(X) + mysvec(V'*[0 prob.obj'/2; prob.obj/2 zeros(n)]*V)'*mysvec(X);



    % obtain the range space
    V = pars.V;
    % V = null(W);
    r = size(V,2);

    X = sdpvar(n+1);
    R = sdpvar(r);

    [~,Asub] = licols([A1 b1]');
    Asub = Asub';
    A1 = Asub(:,1:end-1);
    b1 = Asub(:,end);


    % set up the model
    F = [R>=0; A1*mysvec(X) == b1; A2*mysvec(X) <= b2; X == V*R*V'];

    % the objective function
    obj = mysvec(blkdiag(0,prob.Q))'*mysvec(X) + mysvec([0 prob.obj'/2; prob.obj/2 zeros(n)])'*mysvec(X);



else
    X = sdpvar(n+1);

    [~,Asub] = licols([A1 b1]');
    Asub = Asub';
    A1 = Asub(:,1:end-1);
    b1 = Asub(:,end);


    % remove tiny entries after FR reformulation
    A1(abs(A1)<1e-12) = 0;
    A2(abs(A2)<1e-12) = 0;
    b1(abs(b1)<1e-12) = 0;
    b2(abs(b2)<1e-12) = 0;



    % set up the model
    F = [X>=0; A1*mysvec(X) == b1; A2*mysvec(X) <= b2];

    % the objective function
    obj = mysvec(blkdiag(0,prob.Q))'*mysvec(X) + prob.obj'*X(2:end,1);
end


opts = sdpsettings('solver','mosek','debug',1,'verbose',0,'savesolveroutput',1);
output = optimize(F,obj,opts);


res = output.solveroutput.res;
if pars.FR
    output.r = length(R);
else
    output.r = n + 1;
end
output.lb = value(obj);
output.status = res.rcodestr;
output.solsta = res.sol.itr.solsta;
output.prosta = res.sol.itr.prosta;



if strcmp(output.status,'MSK_RES_OK')
    output.status = 'RES_OK';
elseif strcmp(output.status,'MSK_RES_TRM_STALL')
    output.status = 'STALL';
end

if strcmp(output.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
    output.prosta = 'PD_FEA';
end

% keyboard
% if ~pars.FR
% 
%     params.OutputFlag = 0;
%     % params.TimeLimit = output_SDP{2}.solvertime;
%     output_gurobi = gurobi(prob,params);
% 
%     output.lb
%     output_gurobi.objval
%     keyboard
% end

end


function A1 = compute_VAV(A,V)
[n,r] = size(V);
m = size(A,1);
r1 = nchoosek(r+1,2);
A1 = zeros(m,r1);
for i = 1:m
    temp = V'*mysmat(A(i,:),n)*V;
    temp = (temp+temp')/2;
    A1(i,:) = mysvec(temp);
end

end