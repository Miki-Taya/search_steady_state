%end_timeでの状態の値と定常値の差を返す関数

function plot_generator_state(cnt,tspan,steady_generator_state)

  %パラメータ設定
  y12 = imag(inv(0.085i));  %1-2間送電線のインピーダンス：z12=0.085j
  y23 = imag(inv(0.092i));  %2-3間送電線のインピーダンス：z32=0.092j
  Xd = [1.6;1.4;1.2];
  Xq = [0.25;0.15;0.15];
  B = [y12 -y12 0; -y12 y12+y23 -y23; 0 -y23 y23];  %B：アドミタンス行列Yの虚部であるサセプタンス行列
  Bred = - inv(diag(Xq) - diag(Xq)*B*diag(Xq));


  error = 10.^cnt;
  initial_generator_state = steady_generator_state + [1;1;1;0;0;0;1;1;1]*error;

  delta_star = steady_generator_state(1:3);
  E_star = steady_generator_state(7:9);

  [Pmech_star,Vfield_star] = get_steady_Pmech_Vfield(delta_star,E_star,Bred,Xq,Xd);


  get_dx_nonlinear_ode_wrap = @(t, generator_state) get_dx_nonlinear_ode(t, generator_state, Xd, Xq, Bred, B, Pmech_star, Vfield_star);

  [t_sol generator_state_sol] = ode45(get_dx_nonlinear_ode_wrap, tspan, initial_generator_state);

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
