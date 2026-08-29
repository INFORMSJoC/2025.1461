function output = SDP_FR(prob,pars)
% this functions considers the RLT+PSD relaxation for a given BIP instance
% INPUT: prob is an BIP instance in gurobi format
% OUTPUT: a time needed for computing the exposing vector

% version 1: Shor's SDP relaxation is of the form
% A1(X)=b1
% A2(X)<b2
% where X is psd
%
% form the FR auxiliary problem
% [A1' A2'][u;v] is psd
% v is non-negative
% [b1' b2']y = 0
% trace([A1' A2']y) = 1

% version 2: relaxed Shor's SDP relaxation
% We consider the following SDP problem
% <A1,X> + B1y = b1
% <A2,X> + B2y <= b1

% form the FR auxiliary problem
% [A1' A2']y is psd
% [B1' B2']y is non-negative
% [b1' b2']y = 0
% trace([A1' A2']y) = 1


% choose FR or not
if pars.FR == true
    V = pars.V;
    % there is no effective reduction
    if size(V,1)==size(V,2)
        % then just act as if we do not choose to perform FR
        pars.FR = false;
    end
end

% perform FR based on solver
if strcmp(pars.solvername,'yalmip')
    output = ysolve(prob,pars);
elseif strcmp(pars.solvername,'mosek')
    output = msolve(prob,pars);
end

end

function output = ysolve(prob,pars)

% load the setting
yalmip('clear')

% form the data matrices
[A1,B1,b1,A2,B2,b2] = formSDPdata(prob,pars);
% number of equality and inequalities
m1 = size(A1,1);
m2 = size(A2,1);

% the number of binary and all variables
n = length(prob.obj);
nb = sum(or(prob.vtype=='B',prob.vtype=='I'));
nc = n - nb;

if ismember(pars.sdptype, ["SHOR",'DNN'])
    % when we have a full size matrix
    % set up the problem in yalmip
    u = sdpvar(m1,1);
    v = sdpvar(m2,1);


    if ~pars.FR
        W = sdpvar(n+1);
        F = [W>=0; mysvec(W) == [A1; A2]'*[u;v]; v>=0;u'*b1 + v'*b2 ==0; trace(W) == 1];
    else
        V = pars.V;
        U = null(V');
        r = size(U,2);
        W = sdpvar(r);


        % keyboard
        % OPTION1
        %   % form the equality constraints
        M = [get_UWU(U) -A1' -A2';              % mysvec(xU*W*U') == [A1; A2]'*[u;v]
            zeros(1,nchoosek(r+1,2)) b1' b2';   % u'*b1 + v'*b2 ==0; 
            mysvec(eye(r))' zeros(1,m1+m2)];     % trace(W) == 1
   
        % OPTION2
        % % form the equality constraints
        % M = [get_UWU(U) -A1' -A2';              % mysvec(xU*W*U') == [A1; A2]'*[u;v]
        %     zeros(1,nchoosek(r+1,2)) b1' b2'];   % u'*b1 + v'*b2 ==0; 
            
        % QR is very expensive for larger problems
        % % remove linearly dependent columns
        % [~,M] = licols(M',1e-10);

        % % transpose back
        % M = [M'; ...
        %     mysvec(eye(r))' zeros(1,m1+m2)];     % trace(W) == 1


        M(abs(M)<1e-12) = 0;

         % size(M,1)
        % form the problem
        F = [W>=0; M*[mysvec(W);u;v] == [zeros(size(M,1)-1,1);1]; v>=0;];
        % F = [W>=0; mysvec(U*W*U') == [A1; A2]'*[u;v]; v>=0;u'*b1 + v'*b2 ==0; trace(W) == 1];

        % F = [F;[get_UWU(U) -A1' -A2';zeros(1,nchoosek(r+1,2)) b1' b2';  mysvec(eye(r))' zeros(1,m1+m2)]*[mysvec(W);u;v] == [zeros(size(M,1)-1,1);1]; ];
    end



elseif ismember(pars.sdptype, ["SHOR_B"])
    % the matrix variable is for binary variables only
    % set up the problem in yalmip
    if ~pars.FR
        W = sdpvar(nb+1);
        u = sdpvar(m1,1);
        v = sdpvar(m2,1);
        F = [W>=0; mysvec(W) == [A1; A2]'*[u;v]; [B1;B2]'*[u;v]==0; v>=0;u'*b1 + v'*b2 ==0; trace(W) == 1];
    else
        V = pars.V;
        U = null(V');
        r = size(U,2);
        W = sdpvar(r);
        u = sdpvar(m1,1);
        v = sdpvar(m2,1);
        F = [W>=0; mysvec(U*W*U') == [A1; A2]'*[u;v]; [B1;B2]'*[u;v]==0; v>=0;u'*b1 + v'*b2 ==0; trace(W) == 1];
    end
end

% keyboard
opts = sdpsettings('solver','mosek','debug',1,'verbose',0,'savesolveroutput',1);
% opts = sdpsettings('solver','mosek','debug',1,'verbose',1);
% opts = sdpsettings('solver','mosek','debug',1,'verbose',1,'cachesolvers',1);
% opts.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600;
yalmip_output = optimize(F,[],opts);


% keyboard
res = yalmip_output.solveroutput.res;
output.solvertime = yalmip_output.solvertime;
output.yalmiptime = yalmip_output.yalmiptime;
output.r = size(W,1);

output.status = res.rcodestr;
if strcmp(res.rcodestr,'MSK_RES_ERR_SPACE')
    output.rk = -1;
    output.prosta = 'ERR_SPACE';
    output.solsta = 'ERR_SPACE';
    return
end


% if trace(pars.X*W) > 1e-5
%     output.cs =
% end

% note that yalmip solves the dual of the FR auxiliary problem. thus, if it
% says the dual infeasible, it means the above FR auxiliary problem is
% infeasible
if ismember(res.sol.itr.prosta,['DUAL_INFEASIBLE'])
    % CASE 1: the solver finds certifciate that the problem is infeasible
    res.sol.itr.prosta = 'D_INFEA';
    output.rk = 0;
elseif strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
    % CASE 2: the problem is feasible.
    res.sol.itr.prosta = 'PD_FEA';
    W = value(W);
    if pars.FR
        W = U*W*U';
        W = (W+W')/2;
    end
    output.W = W;
    [Wv,Wd] = eig(W);
    idx = diag(Wd)>1e-4;
    output.rk = sum(idx);
    output.V = Wv(:,~idx);
    
else
    output.rk = -1;
end

% save the information
output.prosta = res.sol.itr.prosta;
output.solsta = res.sol.itr.solsta;

if strcmp(output.status,'MSK_RES_TRM_STALL')
output.status = 'TRM_STALL';
end

if strcmp(output.solsta,'DUAL_INFEASIBLE_CER')
output.solsta  = 'D_INFEA';
end

% keyboard
yalmip('clear')
end




function output = msolve(prob,pars)
solvertime = tic;
[~, res0] = mosekopt('symbcon echo(0)');
symbcon = res0.symbcon;

% form the data matrices
[A1,B1,b1,A2,B2,b2,nb] = formSDPdata(prob,pars);

m1 = length(b1);
m2 = length(b2);
m = m1 +m2;

nc = length(prob.obj) - nb;


% param.MSK_IPAR_LOG = 4;
% param.MSK_SPAR_LOG_FILE_NAME = 'mosek_log.txt';


if pars.FR == false
    r = nb + 1;
    % trace([A1' A2']y) = 1
    n1 = nchoosek(r+1,2);
    M = zeros(r);
    M(tril(true(r))) = 1:n1;
    idx = diag(M);
    M = sparse(1,idx,1,1,n1);
    C = M*[A1' A2'];

    % [B1' B2']y is non-negative
    % [b1' b2']y = 0
    % trace([A1' A2']y) = 1
    probm.a = [B1' B2'; b1' b2';C];
    probm.blc = [zeros(nc+1,1);1];
    probm.buc = [zeros(nc,1); 0;1];

    % [A1' A2']y is psd
    probm.f = sparse([A1' A2']);
    probm.accs = [symbcon.MSK_DOMAIN_SVEC_PSD_CONE nchoosek(r+1,2)];

    probm.blx = [-Inf(m1,1);zeros(m2,1)];
    probm.bux = Inf(m,1);

else
    % R is psd
    % [A1' A2']y - URU'  = 0
    % [B1' B2']y is non-negative
    % [b1' b2']y = 0
    % trace([A1' A2']y) = 1
    V = pars.V;

    U = null(V');
    r = size(U,2);

    % form the mapping URU'
    n1 = nchoosek(nb+2,2);
    r1 = nchoosek(r+1,2);
    URU = zeros(n1,r1);
    k = 0;
    for j = 1:r
        for i = j:r
            k = k + 1;
            if i == j
                R = sparse(i,j,1,r,r);
            else
                R = sparse([i j],[j i],[1 1],r,r);
            end
            temp = U*R*U';
            temp = (temp+temp')/2;
            URU(:,k) = mysvec(temp);
        end
    end

    URUm = cell(n1,1);
    k = 0;
    for j = 1:nb+1
        for i = j:nb+1
            k = k + 1;
            URUm{k} = mysmat(URU(k,:),r);
        end
    end
    % keyboard

    % for the matrix for trace([A1' A2']y) = 1
    n1 = nchoosek(nb+2,2);
    M = zeros(nb+1);
    M(tril(true(nb+1))) = 1:n1;
    idx = diag(M);
    M = sparse(1,idx,1,1,n1);
    C = M*[A1' A2'];

    % URU' - [A1' A2']y = 0
    % [B1' B2']y is non-negative
    % [b1' b2']y = 0
    % trace([A1' A2']y) = 1
    % R >= 0
    probm.bara = get_bara(URUm);
    probm.a = [-A1' -A2';...
        B1' B2';...
        b1' b2';...
        C];
    probm.blc = [zeros(n1,1);zeros(nc,1);0;1];
    probm.buc = [zeros(n1,1);zeros(nc,1);0;1];


    probm.blx = [-Inf(m1,1);zeros(m2,1)];
    probm.bux = [Inf(m1,1); Inf(m2,1)];
end

% solve the problem
try
    mosek_pars = [];
    % mosek_pars.MSK_IPAR_NUM_THREADS = 2;
    [~, res] = mosekopt('minimize info echo(0)', probm,mosek_pars);
    % [~, res] = mosekopt('minimize',probm,mosek_pars);
catch ME
    fprintf('[MOSEK CRASH] MATLAB caught error:\n%s\n', ME.message);
    output.rk = -1;
    output.status = ME.message;
    output.prosta = 'na';
    output.solsta = 'na';
    return
end

% the problem size
output.r = r;

% save the exit flag
output.status = res.rcodestr;
if strcmp(res.rcodestr,'MSK_RES_ERR_SPACE')
    output.rk = -1;
    output.prosta = 'na';
    output.solsta = 'na';
    output.solvertime = toc(solvertime);
    return
end

output.prosta = res.sol.itr.prosta;
output.solsta = res.sol.itr.solsta;


if strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
    % CASE 1: the solver finds certifciate that the problem is infeasible
    output.rk = 0;
elseif strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
    % CASE 2: the problem is feasible.
    if ~pars.FR
        y = res.sol.itr.xx;
        W = mysmat(probm.f*y,nb+1);
        output.rk = sum(eig(W) > 1e-5);
    else
        if ~isempty(res.sol.itr.barx)
            % find out the reduction on psd variable
            R = mysmat(res.sol.itr.barx,r);
            output.rk = sum(eig(R) > 1e-5);
        else
            % no reduction on the psd variable
            output.rk = 0;
        end
    end
else
    output.rk = -1;
    %     res.rcodestr
    %     res.sol.itr.prosta
    % keyboard
end

output.solvertime = toc(solvertime);

% keyboard
% mosek_pars.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-12;
% mosek_pars.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-12;
% mosek_pars.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-12;
% [~, res] = mosekopt('minimize info', probm,mosek_pars);
%
% % dual_psd = ;  % assuming this is the only acc
% Xpsd = mysmat(res.sol.itr.doty,nb+1);

% output.res = res;

% y = res.sol.itr.xx;
% err = sum(min(probm.a*y - probm.blc,0).^2) + sum(max(probm.a*y - probm.buc,0).^2);
% min(eig(mysmat(probm.f*y,nb+1)))
% err = sum(min(y - probm.blx,0).^2) + sum(max(y - probm.bux,0).^2);
%
%
%
% y = res.sol.itr.xx;
% err = 0
% for i = 1:size(probm.a,1)
%     err = err + min(trace(URUm{i}*W)+ probm.a(i,:)*y - probm.blc,0)^2;
%     err = err + max(trace(URUm{i}*W)+ probm.a(i,:)*y - probm.buc,0)^2;
% end

% keyboard
% X = my_testp(prob,pars)


% [y probm.blx probm.bux]`
end





function X = my_testA(prob,pars)


% form the data matrices
[A1,B1,b1,A2,B2,b2] = formSDPdata(prob,pars);

% number of equality and inequalities
m1 = size(A1,1);
m2 = size(A2,1);
nb = sum(or(prob.vtype == 'B',prob.vtype == 'I'));
% set up the problem in yalmip
n = length(prob.obj);
X = sdpvar(n+1);
F = [X>=0; X(1,1) == 1; diag(X(1:nb+1,1:nb+1)) == X(1:nb+1,1); A1*mysvec(X)== b1;A2*mysvec(X) <= b2];


opts = sdpsettings('solver','mosek','debug',1,'verbose',0);
output = optimize(F,[],opts)

X = value(X);



% x = pars.x(:,1)
%
% Aeq*x - beq
% norm(max(Aineq*x - bineq,0))
%
% X0 = [1;x(1:nb)]*[1;x(1:nb)]'
% y0 = x(nb+1:end)
% A1*mysvec(X0) + B1*y0 - b1

end


function X = my_testB(prob,pars)


% form the data matrices
[A1,B1,b1,A2,B2,b2,nb] = formSDPdata(prob,pars);

% number of equality and inequalities
m1 = size(A1,1);
m2 = size(A2,1);

% set up the problem in yalmip
nc = length(prob.obj) - nb;
X = sdpvar(nb+1);
y = sdpvar(nc,1);
F = [X>=0; X(1,1) == 1; diag(X) == X(:,1); A1*mysvec(X) + B1*y == b1;A2*mysvec(X) + B2*y <= b2];


opts = sdpsettings('solver','mosek','debug',1,'verbose',0);
output = optimize(F,[],opts)

X = value(X);



% x = pars.x(:,1)
%
% Aeq*x - beq
% norm(max(Aineq*x - bineq,0))
%
% X0 = [1;x(1:nb)]*[1;x(1:nb)]'
% y0 = x(nb+1:end)
% A1*mysvec(X0) + B1*y0 - b1

end



function M = get_UWU(U)

% get the matrix representation of mysvec(U*W*U')
[n,r] = size(U);
M = zeros(nchoosek(n+1,2),nchoosek(r+1,2));
k = 0;
for j = 1:r
    for i = j:r
        if i == j
            temp = sparse(i,j,1,r,r);
        else
            temp = sparse([i j],[j i],[1/sqrt(2) 1/sqrt(2)],r,r);
        end
        k = k + 1;
        M(:,k) = mysvec(U*temp*U');
    end
end

% find M so that
% svec(U*W*U') = M*svec(W)
% The matrix vector multiplication M*svec(W), means for the column M(:,k),
% we are multiplying svec(W) with 1 in the kth entry,
% this means the matrix W has Wij = 1/sqrt(2)

end
