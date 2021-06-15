function plot_FG_accum_func(steady_generator_state, delta, deltaomega, E, Bred, Xd, Xq, Vfield_star, omega0, M, t_sol, flag_accum_diff)
  
  delta_star = steady_generator_state(1:3);
  deltaomega_star = transpose(steady_generator_state(4:6));
  E_star = steady_generator_state(7:9);
  
  [sol_size,~] = size(t_sol);

  
%F：非線形微分代数方程式系と非線形常微分方程式系は発電機バスのクロン縮約だから、機械サブシステムの方はどちらもこれで同じ
  W_F = zeros(1,sol_size);

  for t = 1:sol_size
    
    W_F(t) = omega0/2 * (deltaomega(t,:)-deltaomega_star) * diag(M) * transpose(deltaomega(t,:)-deltaomega_star);

  end

  W_F = transpose(W_F);



%G：非線形微分代数方程式から発電機バスのクロン縮約をして、非線形常微分方程式系で考えたポテンシャルエネルギー関数(U_G)と蓄積関数(Wred_G)
% nabla は nabla を含む項全体を表す
 
  Ured_G = zeros(1,sol_size);
  Ured_G_star = 0;
  sigma_cos_star = 0;
  nabla = zeros(1,sol_size);

  
  for t = 1:sol_size
      
    for i = 1:3
      sigma_cos = 0;
      sigma_sin = 0;
      
      for j = 1:3
          
          if t == 1
            sigma_cos_star = sigma_cos_star + E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j));
          end
          
        sigma_cos = sigma_cos + E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j));
        sigma_sin = sigma_sin + E(t,j)*Bred(i,j)*sin(delta(t,i)-delta(t,j));
      end
      
          if t == 1 
             Ured_G_star = Ured_G_star + Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E_star(i)*sigma_cos_star; 
          end
         
      Ured_G(t) = Ured_G(t) + Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E(t,i)*sigma_cos;
      nabla(t) = nabla(t) - E(t,i)*sigma_sin*(delta(t,i)-delta_star(i)) + Vfield_star(i)/(Xd(i)-Xq(i))*(E(t,i)-E_star(i));
    end

  end
  
  Ured_G = Ured_G/2;
  Ured_G_star = Ured_G_star/2;
  
  Ured_G = transpose(Ured_G);
  nabla = transpose(nabla);

  
  Wred_G = Ured_G - Ured_G_star - nabla;
  
  figure;
  plot(t_sol, Ured_G)
  title("U^{red}_G")
  
  figure;
  plot(t_sol, W_F)
  title("W_F")
  
  figure;
  plot(t_sol, Wred_G)
  title("W^{red}_G")
  
  figure;
  plot(t_sol,W_F+Wred_G)
  title("W_F + W^{red}_G")

  
  
  % flag_accum_diff == 1 なら [diff( W_F )], [diff( Wred_G )], [diff( W_F + W^{red}_G )] を表示
  if flag_accum_diff == 1
 
    dff_F = diff(W_F);
    [sz,~] = size(dff_F);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_F)
    title("diff( W_F )")
    
    dff_redG = diff(Wred_G);
    [sz,~] = size(dff_redG);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_redG)
    title("diff( Wred_G )")
         
    dff_FredG = diff(W_F+Wred_G);
    [sz,~] = size(dff_FredG);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_FredG)
    title("diff( W_F + W^{red}_G )")
    
  end



end