%end_timeでの状態の値と定常値の差を返す関数

function plot_generator_state(cnt,tspan,initial_generator_state, Pmech_star, Vfield_star)

  %パラメータ設定
  Xq = [0.9360;0.9110;0.6670];
  Xd = [1.5690;1.6510;1.2200];
  Bred = [-6.1331,1.4914,1.6779; 1.4914,-5.9131,2.2693; 1.6779,2.2693,-5.6149];

  error = cnt;
  initial_generator_state = steady_generator_state + [1;1;1;0.0001;0.0001;0.0001;1;1;1]*error

  delta_star = steady_generator_state(1:3);
  E_star = steady_generator_state(7:9);

  [Pmech_star,Vfield_star] = get_steady_Pmech_Vfield(delta_star,E_star,Bred,Xq,Xd)


  get_dx_nonlinear_ode_wrap = @(t, generator_state) get_dx_nonlinear_ode(t, generator_state, Xd, Xq, Bred, Pmech_star, Vfield_star);

  [t_sol, generator_state_sol] = ode45(get_dx_nonlinear_ode_wrap, tspan, initial_generator_state);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);

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
