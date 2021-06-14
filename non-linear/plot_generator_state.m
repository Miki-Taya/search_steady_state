
function t_U = plot_generator_state(error,tspan,steady_generator_state)

  %パラメータ設定
  Xq = [0.9360;0.9110;0.6670];
  Xd = [1.5690;1.6510;1.2200];
  BB = [-6.1331,1.4914,1.6779; 1.4914,-5.9131,2.2693; 1.6779,2.2693,-5.6149];  %BB：アドミタンス行列Yの虚部であるサセプタンス行列
  Bred = - inv(diag(Xq) - diag(Xq)*BB*diag(Xq));

 
  initial_generator_state = steady_generator_state + error

  delta_star = steady_generator_state(1:3);
  E_star = steady_generator_state(7:9);

  [Pmech_star,Vfield_star] = get_steady_Pmech_Vfield(delta_star,E_star,Bred,Xq,Xd);


  get_dx_nonlinear_ode_wrap = @(t, generator_state) get_dx_nonlinear_ode(t, generator_state, Xd, Xq, Bred, Pmech_star, Vfield_star);

  [t_sol, generator_state_sol] = ode45(get_dx_nonlinear_ode_wrap, tspan, initial_generator_state);

  delta = generator_state_sol(:,1:3);
  deltaomega = generator_state_sol(:,4:6);
  E = generator_state_sol(:,7:9);

  [sol_size,~] = size(t_sol);
  
  U = zeros(1,sol_size);
  temp_U = zeros(sol_size,3);
  
  for t = 1:sol_size
  U(t) = 0;
    for i = 1:3
      temp_U(t,i) = 0;
      for j = 1:3
        temp_U(t,i) = temp_U(t,i) + E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j));
      end
      U(t) = U(t) + Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E(t,i)*temp_U(t,i) ;
    end
    U(t) = U(t)/2;
  end
  U = transpose(U);


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
  
  %{
  subplot(4,1,4)
  plot(t_sol, U)
  ylabel('U')
  legend('U')
  %}
  diff_deltaomega = deltaomega(end,:)
  diff_E = E(end,:) - E_star
  
  
  %{
  %蓄積関数:W が半正定関数がどうか調べる
  if all(diff(diff(U)) >= 0)
      disp('W is positive semi-definite.');
  else
      disp('W is not positive semi-definite.');
  end
  %}
  
  
  %--------------------------------------------------------------
  %　必要じゃないなら、ここを消して、関数の返り値や呼び出す側の変数もなくす！
  t_U = [t_sol,U];
  %--------------------------------------------------------------


end