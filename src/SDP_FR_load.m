function SDP_FR_load(idx)
load('mytemp.mat')
output = SDP_FR(prob,parsFR);
clear idx
save('mytemp.mat')


end