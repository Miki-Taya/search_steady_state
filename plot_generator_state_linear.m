%end_timeでの状態の値と定常値の差を返す関数

function plot_generator_state_linear(cnt,tspan,steady_generator_state)

  error = [1;1;1;1;1;1;0;0;0].*(10.^cnt);

  psi = get_psi(steady_generator_state);

  get_dx_linear_ode_wrap = @(tspan, error) get_dx_linear_ode(tspan, error, psi);

  [t_sol, generator_state_sol] = ode45(get_dx_linear_ode_wrap, tspan, error);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,7:9);
  E = generator_state_sol(:,4:6);

  subplot(3,1,1)
  plot(t_sol, delta)
  ylabel('delta')
  legend('delta1','delta2','delta3')

  subplot(3,1,2)
  plot(t_sol, deltaomega)
  ylabel('deltaomega')
  legend('deltaomega1','deltaomega2','deltaomega3')
  
  subplot(3,1,3)
  plot(t_sol, E)
  ylabel('E')
  legend('E1','E2','E3')

end
