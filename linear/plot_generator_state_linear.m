% plotして、end_timeでの状態の値と定常値の差を返す関数

function last = plot_generator_state_linear(error, tspan, psi)

  get_dx_linear_ode_wrap = @(tspan, error) get_dx_linear_ode(tspan, error, psi);

  [t_sol, generator_state_sol] = ode45(get_dx_linear_ode_wrap, tspan, error);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);

  last = generator_state_sol(end,:);
  
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
  
  sgt = sgtitle(['＜線形モデル応答＞　初期値:',num2str(transpose(error))]);
  sgt.FontSize = 10;

end
