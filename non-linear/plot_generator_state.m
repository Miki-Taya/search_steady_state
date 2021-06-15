
function plot_generator_state(error,tspan,steady_generator_state,flag_accum,flag_accum_diff)

  %パラメータ設定
  Xq = [0.9360;0.9110;0.6670];
  Xd = [1.5690;1.6510;1.2200];
  BB = [-6.1331,1.4914,1.6779; 1.4914,-5.9131,2.2693; 1.6779,2.2693,-5.6149];  %BB：アドミタンス行列Yの虚部であるサセプタンス行列
  Bred = - inv(diag(Xq) - diag(Xq)*BB*diag(Xq));
  omega0 = 376.9911;  %系統周波数：60[Htz]*2pi
  M = [100, 18, 12];

 
  initial_generator_state = steady_generator_state + error

  delta_star = steady_generator_state(1:3);
  deltaomega_star = steady_generator_state(4:6);
  E_star = steady_generator_state(7:9);

  [Pmech_star,Vfield_star] = get_steady_Pmech_Vfield(delta_star,E_star,Bred,Xq,Xd);


  get_dx_nonlinear_ode_wrap = @(t, generator_state) get_dx_nonlinear_ode(t, generator_state, Xd, Xq, Bred, Pmech_star, Vfield_star, omega0, M);

  [t_sol, generator_state_sol] = ode45(get_dx_nonlinear_ode_wrap, tspan, initial_generator_state);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);
  
  subplot(1,3,1)
  plot(t_sol, delta)
  yline(delta_star)
  ylabel('delta')
  legend('delta1','delta2','delta3')

  subplot(1,3,2)
  plot(t_sol, deltaomega)
  yline(0)
  ylabel('deltaomega')
  legend('deltaomega1','deltaomega2','deltaomega3')

  subplot(1,3,3)
  plot(t_sol, E)
  yline(E_star)
  ylabel('E')
  legend('E1','E2','E3')
  %axis([0 100 3.21 3.24]) %軸の範囲指定　x軸[0 100]  y軸[3.21 3.24]
  sgtitle(['初期誤差: ',num2str(transpose(error))]);
  
  %最終値と定常値との差を表示したいなら ; 外して
  diff_deltaomega = deltaomega(end,:) - deltaomega_star;
  diff_E = E(end,:) - E_star;
  
  

  if flag_accum == 1
      plot_FG_accum_func(steady_generator_state, delta, deltaomega, E, Bred, Xd, Xq, Vfield_star, omega0, M, t_sol, flag_accum_diff)
  end

  
  
  
  %--------------------------------------------------------------
  %　必要じゃないなら、ここを消して、関数の返り値や呼び出す側の変数もなくす！
  %t_U = [t_sol,U];
  %--------------------------------------------------------------


end
