% Run every table-producing experiment using its complete problem list.
% This requires all MIPLIB files referenced by prob_list_A.mat and
% prob_list_B.mat. Use smoke_test for a partial run with the included data.
for EXPnr = 2:5
    test_primalFR(EXPnr)
end
