function phiAltVEFa=solveAdjointAltVEF
global dat results
[phiAltVEFa]=solve_VEF(~dat.forward_flux,results.Ea,results.Ebda);
return
end