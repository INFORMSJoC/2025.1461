% Run the first few instances of every table-producing experiment.
% The included 48-instance data subset contains all files needed here.
% Experiment 5 can still take at least 600 seconds on an instance because
% the instance limit does not change the per-instance solver time limits.

instance_limit = 3;
for EXPnr = 2:5
    test_primalFR(EXPnr, instance_limit)
end
