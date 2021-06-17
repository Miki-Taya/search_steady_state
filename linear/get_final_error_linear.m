%end_timeでの状態の値と定常値の差を返す関数

function final_error = get_final_error_linear(cnt,tspan, psi)

  error = [1;1;1;0;0;0;1;1;1]*cnt;

  get_dx_linear_ode_wrap = @(tspan, error) get_dx_linear_ode(tspan, error, psi);

  [~, generator_state_sol] = ode45(get_dx_linear_ode_wrap, tspan, error);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);

  final_error = transpose([delta(end,:); deltaomega(end,:); E(end,:)]);

end
