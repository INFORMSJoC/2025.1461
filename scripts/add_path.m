function folder_path = add_path()
% Add the repository and locally installed dependencies to the MATLAB path.

script_path = fileparts(mfilename('fullpath'));
folder_path = [fileparts(script_path) filesep];
addpath(genpath(folder_path));

% Optional environment variables for third-party MATLAB packages.
% Each variable should contain the directory that is to be added to the
% MATLAB path. See the repository README for configuration examples.
dependency_variables = {
    'YALMIP_ROOT'
    'MOSEK_MATLAB_ROOT'
    'GUROBI_MATLAB_ROOT'
};

for i = 1:numel(dependency_variables)
    variable_name = dependency_variables{i};
    dependency_path = getenv(variable_name);

    if isempty(dependency_path)
        continue
    end

    if isfolder(dependency_path)
        addpath(genpath(dependency_path));
    else
        warning('add_path:InvalidDependencyPath', ...
            '%s does not identify an existing directory: %s', ...
            variable_name, dependency_path);
    end
end
end
