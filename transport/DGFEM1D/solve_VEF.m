function phi = solve_VEF(forward,E,Ebd)

%[A,rhs] = build_VEF_system(forward,E,Ebd);
[A,rhs] = build_VEF_system_alt(forward,E,Ebd);

phi = A\rhs;

end