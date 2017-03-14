function phi = solve_VEF_math_adjoint(forward,E,Ebd)

if forward
    error('in %s, we have be in adjoint mode',mfilename);
end

[A,rhs] = build_VEF_system_alt(forward,E,Ebd);
A=A';

phi = A\rhs;

end