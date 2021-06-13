function dx = get_dx_linear_ode(t, x, psi)

  %定義：dx/dt = [ ddelta/dt;  ddeltaomega/dt; dE/dt]
  dx = psi*x;

end