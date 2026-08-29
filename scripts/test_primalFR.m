function test_primalFR(EXPnr,instance_limit)
%TEST_PRIMALFR Run a facial-reduction experiment.
%
%   test_primalFR(EXPnr) runs the configured number of instances.
%   test_primalFR(EXPnr,instance_limit) runs at most the first
%   instance_limit instances. For example, test_primalFR(2,3) runs the
%   first three instances of the DNN experiment.
%   Result MAT files and diary logs are saved in the repository's results
%   subfolder.
%
% This file tests affine generation.
% profile on
% initialization
narginchk(1,2);
validateattributes(EXPnr, {'numeric'}, ...
    {'scalar','integer','>=',1,'<=',6}, mfilename, 'EXPnr', 1);
if nargin < 2 || isempty(instance_limit)
    instance_limit = Inf;
else
    validateattributes(instance_limit, {'numeric'}, ...
        {'scalar','integer','positive'}, mfilename, 'instance_limit', 2);
end

script_path = fileparts(mfilename('fullpath'));
restoredefaultpath;
addpath(script_path);
diary off
folder_path = add_path(); % Add the current folder and its all subfolders to the path.

% experiment number
% EXPnr = 1; % comprehensive test of different primal FR methods
% EXPnr = 2; % Applying a specific primal FR method followed by FR, DNN
% EXPnr = 3; % Applying a specific primal FR method followed by FR, SHOR
% EXPnr = 4; % Applying a specific primal FR method followed by FR, SHOR V
% EXPnr = 5; % Applying a specific primal FR method followed by FR, and SDP

% record setting
record_result = true;
% record_result = false;

% new session
% new_session = true;
new_session = false;

% global settings
DEBUG = true;
DEBUG = false;

% get the title and keep generated files out of the source directories
results_path = fullfile(folder_path, 'results');
if record_result && ~isfolder(results_path)
    mkdir(results_path);
end
fname = fullfile(results_path, get_title(EXPnr));
if record_result
    diary([fname '.txt']);
end

fprintf('Test at %s\n',string(datetime))

% list of primal FR implementation methods
method_names = {'OPT','FEA','OPT(R)','FEA(R)','OPT(W)','FEA(W)','OPT(RW)','FEA(RW)'};
% OPT = optimization based
% FEA = feasibility based
% R = randomized obj
% W = warm_start generating 0-1 solutions


% choose the methods
if EXPnr == 1
    % pick some methods, and perform feasible solution generation
    methods = 1:length(method_names);
    % methods = 7;
    var_mode = 'A';        % variable mode
    mmax = 106;
elseif EXPnr == 2
    methods = 7;           % pick one method, perform FRA
    var_mode = 'A';        % variable mode
    solvername = 'yalmip'; % choose solver for FR
    sdptype = 'DNN';
    mmax = 26;
elseif EXPnr == 3
    methods = 7;           % pick one method, perform FRA
    var_mode = 'A';        % variable mode
    solvername = 'yalmip'; % choose solver for FR
    sdptype = 'SHOR';
    mmax = 50;
elseif EXPnr == 4
    methods = 7;           % pick one method, perform FRA
    var_mode = 'B';        % variable mode
    solvername = 'yalmip'; % choose solver for FR
    sdptype = 'SHOR_B';
    mmax = 106;
elseif EXPnr == 5
    methods = 7;           % pick one method, perform FRA
    var_mode = 'A';        % variable mode
    solvername = 'yalmip'; % choose solver for FR
    % sdptype = 'DNN';
    % mmax = 16;

    sdptype = 'SHOR';
    mmax = 100;
elseif EXPnr == 6
    methods = 7;           % pick one method, perform FRA
    var_mode = 'A';        % variable mode
    solvername = 'yalmip'; % choose solver for FR
    sdptype = 'SHOR';
    mmax = 100;
end



% load the problem list
if strcmp(var_mode,'A')
    load('prob_list_A.mat','prob_list');
elseif strcmp(var_mode,'B')
    load('prob_list_B.mat','prob_list');
end


% start the test
np = length(prob_list); % number of problems
mmax = min([mmax,instance_limit,np]);
fprintf('Experiment %d will run the first %d instance(s).\n',EXPnr,mmax)

m = length(methods);    % number of methods
output = cell(1,m);
output_FR = cell(2,1);
output_table = cell(np,m+1);
output_SDP = cell(2,1);
output_gurobi = [];

% if EXPnr == 3
% load result_EXP3_20250627_001816a.mat
% mmin = 37;
% else
%     mmin = 1;
% end

% load result_EXP2_20250619_161157a.mat
% load result_EXP4_20250621_163043a.mat
% load result_EXP5_20250625_191632a.mat
% print the title
for i = 1:mmax
    % print the time
    fprintf('Current time is %s\n',string(datetime))

    % load the problem and make some corrections
    prob = gurobi_read([folder_path 'data/' prob_list{i}.name]);

    % make simple adjustments
    prob = prob_correction(prob,prob_list{i});

    % print the title and the problem summary
    print_title(methods,method_names)
    print_prob(i,prob,prob_list)

    % generate feasible solutions using different methods
    for j = 1:m
        % get the parameters
        pars = get_pars(method_names{j},prob,DEBUG,var_mode);
        % find the feasible solutions
        output{j} = get_feasible_sols(prob,pars);
    end


    % if we solve the FR auxiliary problem
    if EXPnr >= 2
        % we choose only one method to perform FR
        if length(methods) > 1
            error('pls choose only one method when testing FR\n')
        end

        output_FR{1} = FRfields(0);
        output_FR{2} = FRfields(0);
        output_gurobi = [];

        % if we find at least two feasible solution, then perform FRAs
        if output{1}.refine.size > 1
            % obtain the feasible solutions
            V = output{1}.refine.V;

            if EXPnr <= 4
                % standard FRA
                output_FR{1} = standardFRA(prob,solvername,sdptype,new_session);
            end

            % primal FRA
            output_FR{2} = primalFRA(prob,V,output{1}.aff,pars,solvername,sdptype,new_session);

            % solve SDP
            if EXPnr >= 5
                output_SDP{1} = SDPfields(0);
                output_SDP{2} = SDPfields(0);
                % if there is a considerable reduction in the problem
                if and(output_FR{2}.rk > 0,output_FR{2}.rk/length(prob.obj)>.1)


                    % add a quadratic objective
                    probq = add_quad_obj(prob);

                    % solve the original formulation
                    parsSDP = [];
                    parsSDP.FR = false;
                    parsSDP.sdptype = sdptype;
                    output_SDP{1} = solve_SDP(probq,parsSDP);

                    % solve the FR formulation
                    parsSDP.FR = true;
                    parsSDP.V = output_FR{2}.V;
                    parsSDP.sdptype = sdptype;
                    output_SDP{2} = solve_SDP(probq,parsSDP);

                    % solve by gurobi
                    params.OutputFlag = 0;
                    params.TimeLimit = max(output_SDP{2}.solvertime,600);
                    output_gurobi = gurobi(probq,params);
                end
            end
        end
    end

    % save information for latex table
    output_table(i,:) = save_results(output,output_FR,output_SDP,output_gurobi,prob,methods,EXPnr);


    % print the result
    print_result(output_table(i,:),methods,EXPnr);

    % save the results
    if record_result
        save([fname '.mat'])
    end

    % stop solving more instances
    if EXPnr == 2
        if output_FR{1}.solvertime > 3600
            break
        end
    end

    % stop solving more instances
    if EXPnr >= 5
        if output_SDP{1}.solvertime > 3600
            break
        end
    end
end

fprintf('The experiment is finished at %s\n',string(datetime))
fprintf('\n\n\n\n\n\n')
diary off

% profile off
% profile viewer
% keyboard
end



function output_SDP = SDPfields(k)
% initialize the output
if k == 0
    output_SDP.solvertime = -1;
    output_SDP.yalmiptime = -1;
    output_SDP.r = -1;
    output_SDP.lb = -Inf;
    output_SDP.status = 'na';
    output_SDP.prosta = 'na';
    output_SDP.solsta = 'na';
elseif k == 1
    % if there is no reduction,
    % then slater's condition is proved without FRA
    output_SDP.r = 0;
    output_SDP.solvertime = 0;
    output_SDP.status = 'proved';
end


end




function output_FR = FRfields(k)
% initialize the output
if k == 0
    output_FR.solvertime = -1;
    output_FR.yalmiptime = -1;
    output_FR.rk = -1;
    output_FR.r = -1;
    output_FR.status = 'na';
    output_FR.prosta = 'na';
    output_FR.solsta = 'na';
elseif k == 1
    % if there is no reduction,
    % then slater's condition is proved without FRA
    output_FR.rk = 0;
    output_FR.r = 0;
    output_FR.solvertime = 0;
    output_FR.status = 'proved';
end


end

function output_FR = standardFRA(prob,solvername,sdptype,new_session)

output = FRfields(0);

% perform standard FRA
parsFR = [];
parsFR.FR = false;
parsFR.solvername = solvername;
parsFR.sdptype = sdptype;
if new_session
    % new session to avoid memory issues
    save('mytemp.mat')
    matlab_executable = fullfile(matlabroot, 'bin', 'matlab');
    cmd = sprintf('"%s" -nojvm -batch "SDP_FR_load(1)"', ...
        matlab_executable);
    system(cmd);
    load('mytemp.mat')
    output_FR = output;
else
    % run in the same matlab session, may encounter memory issues.
    output_FR = SDP_FR(prob,parsFR);
end

end


function output_FR = primalFRA(prob,V,output,pars,solvername,sdptype,new_session)

% perform the primal FRA
V1 = get_face(V,pars);% get the face
if size(V1,1) == size(V1,2)
    % Case 1: the feasible solutions are full-dimenional
    % then slater's condition is proved without FRA
    output_FR = FRfields(1);
    output_FR.V = eye(size(V1,2));
elseif strcmp(output.status,'SUCC_SIZE')
    % Case 2: the feasible solutions match the dimension of affine set
    % upper bound, and thus slater's condition is proved without FRA
    output_FR = FRfields(1);
    % the feasible solutions span V
    output_FR.V = V1;
    output_FR.rk = size(V1,1)-size(V1,2);
else
    % Case 3: we do not know if these feasible solutions span affine hull
    % run the facial reduction algorithm
    output = FRfields(0);
    % parameters
    parsFR = [];
    parsFR.FR = true;
    parsFR.V = V1;
    parsFR.solvername = solvername;
    parsFR.sdptype = sdptype;

    % solve the primal FRA problem
    if new_session
        % new session to avoid memory issues
        save('mytemp.mat')
        matlab_executable = fullfile(matlabroot, 'bin', 'matlab');
        cmd = sprintf('"%s" -nojvm -batch "SDP_FR_load(2)"', ...
            matlab_executable);
        system(cmd);
        load('mytemp.mat')
        output_FR = output;
    else
        % run in the same matlab session, may encounter memory issues.
        output_FR = SDP_FR(prob,parsFR);
    end
end

end

function output_table = save_results(output,output_FR,output_SDP,output_gurobi,prob,methods,EXPnr)
% the output of the algorithm is saved in output_table
% output_table.fea := affinely independent solutions
% output_table.refine := feasible solutions after refinement
% output_table.FR1 := standard FRA
% output_table.FR2 := primal FRA

m = length(methods);
output_table = cell(1,m+1);

% save problem information
output_table{end} = prob_info(prob);

% save the experiment results
for j = 1:m
    output_table{j}.implicit = output{j}.implicit;
    output_table{j}.warmup = output{j}.warmup;
    output_table{j}.aff = output{j}.aff;
    output_table{j}.refine = output{j}.refine;

    if EXPnr >= 2
        output_table{j}.FR1 = output_FR{1};
        output_table{j}.FR2 = output_FR{2};
        output_table{j}.SDP1 = output_SDP{1};
        output_table{j}.SDP2 = output_SDP{2};
        output_table{j}.GU = output_gurobi;
    end
end


% remove some big fields to save space
for j = 1:m
    if ~isempty(output_table{j})
        if isfield(output_table{j}.warmup,'V')
            output_table{j}.warmup = rmfield(output_table{j}.warmup,'V');
        end

        if isfield(output_table{j}.refine,'V')
            output_table{j}.refine = rmfield(output_table{j}.refine,'V');
        end

        if isfield(output_table{j}.aff,'V')
            output_table{j}.aff = rmfield(output_table{j}.aff,'V');
        end

        % if EXPnr >= 2
        %     if isfield(output_table{j}.FR1,'W')
        %         output_table{j}.FR1 = rmfield(output_table{j}.FR1,'W');
        %     end
        %     if isfield(output_table{j}.FR2,'W')
        %         output_table{j}.FR2 = rmfield(output_table{j}.FR2,'W');
        %     end
        % end
    end
end

end


function print_result(output,methods,EXPnr)
% print the results
methods = logical(methods);

% each line of the table 'name','format','field','field'
contents = {'Size','d','aff','size';...
    'Rank','d','aff','rk';...
    'Size(max)','d','implicit','K';...
    'warmup','d','warmup','size';...
    '#Sub','d','aff','nr_sub';...
    'Time','.2f','aff','time';...
    'Status','s','aff','status';...
    'cond','.2e','aff','cond';...
    'blank','','','';...
    'Size','d','refine','size';...
    'Rank','d','refine','rk';...
    '#Sub','d','refine','nr_sub';...
    'Time','.2f','refine','time';...
    'cond','.2e','refine','cond';};

if EXPnr >= 2
    contents = [contents;...
        {'blank','','','';...
        'Size','d','FR1','r';...
        'Rank','d','FR1','rk';...
        'Time','.2f','FR1','solvertime';...
        'Status','s','FR1','status';...
        'prosta','s','FR1','prosta';...
        'solsta','s','FR1','solsta';...
        'blank','','','';...
        'Size','d','FR2','r';...
        'Rank','d','FR2','rk';...
        'Time','.2f','FR2','solvertime';...
        'Ytime','.2f','FR2','yalmiptime';...
        'Status','s','FR2','status';...
        'prosta','s','FR2','prosta';...
        'solsta','s','FR2','solsta'}];
end


if EXPnr >= 5
    contents = [contents;...
        {'blank','','','';...
        'Size','d','SDP1','r';...
        'L.B.','.2f','SDP1','lb';...
        'Time','.2f','SDP1','solvertime';...
        'Ytime','.2f','SDP1','yalmiptime';...
        'Status','s','SDP1','status';...
        'prosta','s','SDP1','prosta';...
        'solsta','s','SDP1','solsta';...
        'blank','','','';...
        'Size','d','SDP2','r';...
        'LB','.2f','SDP2','lb';...
        'Time','.2f','SDP2','solvertime';...
        'Ytime','.2f','SDP2','yalmiptime';...
        'Status','s','SDP2','status';...
        'prosta','s','SDP2','prosta';...
        'solsta','s','SDP2','solsta';...
        'blank','','','';...
        'LB','.2f','GU','objval';...
        'Time','.2f','GU','runtime';...
        'Status','s','GU','status'}];
end


for l = 1:size(contents,1)

    % left vertical line
    if l > 1
        fprintf('%-90s%-5s','','|')
    end

    % if it is a blank line
    if strcmp(contents{l,1},'blank')
        fprintf('%s|\n',char(ones(1,10+length(methods)*10)*'-'))
        continue
    end

    % label
    fprintf('%-10s',contents{l,1})

    for ii = 1:length(methods)
        if ~isfield(output{ii},contents{l,3})
            fprintf('%-10s',' ')
            continue
        end
        if ~isfield(output{ii}.(contents{l,3}),contents{l,4})
            fprintf('%-10s',' ')
            continue
        end
        fprintf(['%-10' contents{l,2}],output{ii}.(contents{l,3}).(contents{l,4}))
    end


    % right vertical line
    fprintf('%-5s\n','|')
end

end

function print_title(methods,method_names)
% print the title
fprintf('%s\n',char(ones(1,106 + length(methods)*10)*'-'))
fprintf('%-5s%-25s%-10s%-10s%-10s%-10s%-10s%-10s','Nr','Name','Level','C','Cons','Var','Bvar','Size')
fprintf('%-5s','|')
fprintf('%-10s','')
for i = 1:length(methods)
    fprintf('%-10s',method_names{methods(i)})
end
fprintf('|')

fprintf('\n')

end

function print_prob(i,prob,prob_list)
C = norm([prob.A(:);prob.rhs;prob.lb(isfinite(prob.lb));prob.ub(isfinite(prob.ub))],'fro');
% find a maximal affinely independent vectors
fprintf('%-5d',i)
fprintf('%-25s',prob.modelname)
fprintf('%-10s',prob_list{i}.level)
fprintf('%-10.1e',C)
fprintf('%-10d',size(prob.A,1))
fprintf('%-10d',length(prob.obj))
fprintf('%-10d',sum(get_bvar(prob)))
fprintf('%-10d',prob_list{i}.psize)
fprintf('%-5s','|')
end



function pinfo = prob_info(prob)
% save information about the problem instances for printing
pinfo.name = prob.modelname;        % model name
pinfo.n = length(prob.obj);         % number of variables
pinfo.nb = sum(or(prob.vtype=='I',prob.vtype=='B')); % binary variables
pinfo.m = length(prob.sense);       % total constraints
pinfo.m1 = sum(prob.sense=='=');    % equality constraints
pinfo.m2 = sum(prob.sense~='=');    % inequality constraints
pinfo.nf = sum(prob.lb == prob.ub); % number of fixed variables
pinfo.nl = sum(isfinite(prob.lb)); % number of lower bounds
pinfo.nu = sum(isfinite(prob.ub)); % number of lower bounds
pinfo.vtype = prob.vtype;    % variable type, continuous/binary/integer
end


function pars = get_pars(method_name,prob,DEBUG,var_mode)
% choose a method for generating feasible solutions
% '#OPT','$FEA','#OPT(R)','#FEA(R)','#OPT(W)','$FEA(W)','#OPT(RW)','#FEA(RW)','#OPT(WS)','#OPT(RWS)','#OPT(RWT)','#FEA(RWT)','#OPT(RST)'

% global settings
% var_mode = 'A'; % all variables
% var_mode = 'B'; % binary variables
refine = true;

% get bivnary variable indices
bvar = get_bvar(prob);
n = length(prob.obj);
if strcmp(var_mode,'A')
    max_time = n;
elseif strcmp(var_mode,'B')
    max_time = sum(bvar);
end
sub_time = 5;

% keyboard
pars = [];
switch method_name
    case 'OPT'
        pars.DEBUG = DEBUG;
        pars.mode = 'opt';
        pars.var_mode = var_mode;
        pars.rand = false;
        pars.warm_up = false;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;
    case 'FEA'

        pars.DEBUG = DEBUG;
        pars.mode = 'fea';
        pars.var_mode = var_mode;
        pars.rand = false;
        pars.warm_up = false;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'OPT(R)'

        pars.DEBUG = DEBUG;
        pars.mode = 'opt';
        pars.var_mode = var_mode;
        pars.rand = true;
        pars.warm_up = false;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'FEA(R)'

        pars.DEBUG = DEBUG;
        pars.mode = 'fea';
        pars.var_mode = var_mode;
        pars.rand = true;
        pars.warm_up = false;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'OPT(W)'

        pars.DEBUG = DEBUG;
        pars.mode = 'opt';
        pars.var_mode = var_mode;
        pars.rand = false;
        pars.warm_up = true;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'FEA(W)'

        pars.DEBUG = DEBUG;
        pars.mode = 'fea';
        pars.var_mode = var_mode;
        pars.rand = false;
        pars.warm_up = true;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'OPT(RW)'

        pars.DEBUG = DEBUG;
        pars.mode = 'opt';
        pars.var_mode = var_mode;
        pars.rand = true;
        pars.warm_up = true;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;

    case 'FEA(RW)'

        pars.DEBUG = DEBUG;
        pars.mode = 'fea';
        pars.var_mode = var_mode;
        pars.rand = true;
        pars.warm_up = true;
        pars.solpool = false;
        pars.max_time = max_time;
        pars.sub_time = sub_time;
        pars.refine = refine;
    otherwise
        fprintf('the method %s does not exist.\n',method_name)
        error('exit now\n')
end


if strcmp(pars.var_mode,'B')
    pars.bvar = bvar;
end

end

function fname = get_title(EXPnr)

title = char(datetime('now','format','yyyy-MM-dd HH:mm:ss'));
title(title=='-') = [];
title(title==' ') = '_';
title(title==':') = [];

fname = ['result_EXP' num2str(EXPnr) '_' title];

end

function prob = add_quad_obj(prob)
% % generate random cost matrix
n = length(prob.obj);
Q = abs(randn(n));
Q(rand(n)<1/4) = 0;
Q = (tril(Q)+tril(Q,-1)')/2;
Q = Q + eye(n)*10;
nb = sum(or(prob.vtype=='I',prob.vtype=='B')); % binary variables
c = abs(randn(n-nb,1));
prob.Q = sparse(Q);
prob.obj = [zeros(nb,1);c];



end
