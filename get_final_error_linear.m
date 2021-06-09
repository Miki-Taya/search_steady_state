%end_timeでの状態の値と定常値の差を返す関数

function final_error = get_final_error_linear(cnt,tspan,steady_generator_state)

  error = [1;-1;1;;0;0;0;1;1;1].*(10.^cnt);

  psi = get_psi(steady_generator_state);

  get_dx_linear_ode_wrap = @(tspan, generator_state) get_dx_linear_ode(tspan, generator_state, psi);

  [t_sol, generator_state_sol] = ode45(get_dx_linear_ode_wrap, tspan, error);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);

  final_error = transpose([delta(end,:); deltaomega(end,:); E(end,:)]);

end
