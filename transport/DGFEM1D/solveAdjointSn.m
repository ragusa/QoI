function [phia,Ea,Ebda,psia]=solveAdjointSn
global dat IO_opts 
[phia,Ea,Ebda,psia]=solve_transport(dat.adjoint_flux,dat.do_dsa,IO_opts.console_io);
return
end