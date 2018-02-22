function phi = solve_VEF_math_adjoint(forward,E,Ebd)
    % adjoint mode
    % passing adjoint keyword to build adjoint src and use adjoint BCs
    [A,rhs] = build_VEF_system_alt(forward,E,Ebd);
    % the matrix is transposed to get the adjoint
    A=A';
    % solve for the adjoint flux
    phi = A\rhs;

end