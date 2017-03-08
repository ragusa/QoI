function phi = solve_VEF(forward,E)

[A,rhs] = build_VEF_system(forward,E);

phi = A\rhs;

end