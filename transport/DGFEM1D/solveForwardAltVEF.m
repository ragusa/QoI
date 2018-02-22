function phiAltVEF=solveForwardAltVEF
global npar dat snq IO_opts results
[phiAltVEF]=solve_VEF_math_adjoint(dat.forward_flux,results.Ea,results.Ebda);
return
end