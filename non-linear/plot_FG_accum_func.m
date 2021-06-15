function plot_FG_accum_func(steady_generator_state, delta, deltaomega, E, Bred, Xd, Xq, Vfield_star, omega0, M, t_sol, flag_accum_diff)
  
  delta_star = steady_generator_state(1:3);
  deltaomega_star = transpose(steady_generator_state(4:6));
  E_star = steady_generator_state(7:9);
  
  [sol_size,~] = size(t_sol);

  
%F
  W_F = zeros(1,sol_size);

  for t = 1:sol_size
    
    W_F(t) = omega0/2 * (deltaomega(t,:)-deltaomega_star) * diag(M) * transpose(deltaomega(t,:)-deltaomega_star);

  end

  W_F = transpose(W_F);



%G nabla は nabla を含む項全体を表す
 
  Ured = zeros(1,sol_size);
  Ured_star = 0;
  sigma_cos = zeros(1,3);
  sigma_sin = zeros(1,3);
  temp_Ured_star = 0;
  nabla = zeros(1,sol_size);
  
  for t = 1:sol_size
      
    for i = 1:3
      
      for j = 1:3
          
          if t == 1
            temp_Ured_star = temp_Ured_star + E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j));
          end
          
        sigma_cos(i) = sigma_cos(i) + E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j));
        sigma_sin(i) = sigma_sin(i) + E(t,j)*Bred(i,j)*sin(delta(t,i)-delta(t,j));
      end
      
          if t == 1
             Ured_star = Ured_star + Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E_star(i)*temp_Ured_star;
          end
            
      Ured(t) = Ured(t) + Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E(t,i)*sigma_cos(i);
      nabla(t) = nabla(t) - E(t,i)*sigma_sin(i)*(delta(t,i)-delta_star(i)) + Vfield_star(i)/(Xd(i)-Xq(i))*(E(t,i)-E_star(i));
    end
    
    sigma_cos = zeros(1,3);
    sigma_sin = zeros(1,3);
  end
  
  Ured(t) = Ured(t)/2;
  Ured_star = Ured_star/2;
  
  Ured = transpose(Ured);
  nabla = transpose(nabla);

  
  Wred_G = Ured - Ured_star - nabla;
  
  figure;
  plot(t_sol, Ured)
  title("Ured")
  
  figure;
  plot(t_sol, W_F)
  title("W_F")
  
  figure;
  plot(t_sol, Wred_G)
  title("Wred_G")
  
  
  if flag_accum_diff == 1
    figure;
    plot(t_sol,W_F+Wred_G)
    title("W_F + Wred_G")
  end



end