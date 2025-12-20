function f = poisson_RHS(S,rhs,LHS_screen)
f = 4 * pi * (rhs)/S.epsilon_elec;
f = LHS_screen*f;
end