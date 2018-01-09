function [phi,E,Ebd,psi]=solveForwardSn
global dat IO_opts 
[phi,E,Ebd,psi]=solve_transport(dat.forward_flux,dat.do_dsa,IO_opts.console_io);
return
end