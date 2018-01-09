function phiVEFa=solveAdjointVEF
global npar dat snq IO_opts results
[phiVEFa]=solve_VEF_math_adjoint(dat.adjoint_flux,results.E,results.Ebd);
return
end