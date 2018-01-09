function phiVEF=solveForwardVEF
global dat results
[phiVEF]=solve_VEF(dat.forward_flux,results.E,results.Ebd);
return
end