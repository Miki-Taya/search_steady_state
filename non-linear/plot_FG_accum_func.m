function plot_FG_accum_func(tspan,steady_generator_state, delta, deltaomega, E, Bred, Xd, Xq, Vfield_star, omega0, M, t_sol, flag_accum_diff)
  
  delta_star = steady_generator_state(1:3);
  deltaomega_star = transpose(steady_generator_state(4:6));
  E_star = steady_generator_state(7:9);
  
  [sol_size,~] = size(t_sol);

  
%F：非線形微分代数方程式系と非線形常微分方程式系は発電機バスのクロン縮約だから、機械サブシステムの方はどちらもこれで同じ
  W_F = zeros(sol_size,1);

  for t = 1:sol_size
    
    W_F(t) = omega0/2 * (deltaomega(t,:)-deltaomega_star) * diag(M) * transpose(deltaomega(t,:)-deltaomega_star);

  end

%平衡点で 0 になる確認  
  W_F_star = zeros(sol_size,1);

  for t = 1:sol_size
    
    W_F_star(t) = omega0/2 * (deltaomega_star-deltaomega_star) * diag(M) * transpose(deltaomega_star-deltaomega_star);

  end  
  
  
  


%G：非線形微分代数方程式から発電機バスのクロン縮約をして、非線形常微分方程式系で考えたポテンシャルエネルギー関数(U_G)と蓄積関数(Wred_G)

 

% 1.1 matrix を for文で表した　こっちはよくない


  % nabla は nabla を含む項全体を表す
  %{
  
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
  
  %}
  

  
% 1.2 平衡点を代入したら W は 0 にならなければならない。そのテスト。for ver.

  %{
  
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
          
        sigma_cos = sigma_cos + E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j));
        sigma_sin = sigma_sin + E_star(j)*Bred(i,j)*sin(delta_star(i)-delta_star(j));
      end
      
          if t == 1 
             Ured_G_star = Ured_G_star + Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E_star(i)*sigma_cos_star; 
          end
         
      Ured_G(t) = Ured_G(t) + Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E_star(i)*sigma_cos;
      nabla(t) = nabla(t) - E_star(i)*sigma_sin*(delta_star(i)-delta_star(i)) + Vfield_star(i)/(Xd(i)-Xq(i))*(E_star(i)-E_star(i));
    end

  end
  
  Ured_G = Ured_G/2;
  Ured_G_star = Ured_G_star/2;
  
  Ured_G = transpose(Ured_G);
  nabla = transpose(nabla);

  
  Wred_G = Ured_G - Ured_G_star - nabla;
  %}  

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
  
  
% 2.1 matrix は matrix のまま計算... nablaU * x の計算
  
  Ured_G = zeros(sol_size,1);
  Ured_G_star = 0;
  trans_nablaU = zeros(1,6);
  x_G = zeros(6,1);
  Wred_G = zeros(sol_size,1);


  
  for t = 1:sol_size
      
        for i = 1:3
              %電気サブシステムGの出力: y_G = -Ei * Σ(j=1,2,3) Ej * Bredij * sinδij
              y_G = 0;
              
              % i 番目の Ured_G を計算する.　　Ured_Gi_star は平衡点での値
              Ured_Gi = (Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i)))) /2;
              Ured_Gi_star = (Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i)))) /2;
              
              for j = 1:3
                  
                          % star の値は t によらないため、t == 1 のときだけ計算
                          if t == 1
                            Ured_Gi_star = Ured_Gi_star + (E_star(i)*E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j))) /2;
                          end
                          
                  Ured_Gi = Ured_Gi + (E(t,i)*E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j))) /2;
                  y_G = y_G - E(t,i) * E(t,j)*Bred(i,j)*sin(delta(t,i)-delta(t,j));
              
              end
              
                          if t == 1 
                             Ured_G_star = Ured_G_star + Ured_Gi_star; 
                          end
                          
              % Ured_G = Ured_G1 + Ured_G2 + UredG3
              Ured_G(t) = Ured_G(t) + Ured_Gi;
              
              % trans_nablaU = [y_G1,y_G2,yG3,Vfield_star(1)/(Xd(1)-Xq(1)),Vfield_star(i)/(Xd(2)-Xq(2)),Vfield_star(3)/(Xd(3)-Xq(3))]
              trans_nablaU(i) = y_G;
              trans_nablaU(i+3) = Vfield_star(i)/(Xd(i)-Xq(i));
              
              % x_G = [delta(t,1)-delta_star(1); delta(t,2)-delta_star(2); delta(t,3)-delta_star(3); E(t,1)-E_star(1); E(t,2)-E_star(2); E(t,3)-E_star(3)];
              x_G(i) = delta(t,i)-delta_star(i);
              x_G(i+3) = E(t,i)-E_star(i);

        end
        Wred_G(t) = Ured_G(t) - Ured_G_star - trans_nablaU * x_G;
        if mod(t,100) == 0
            Ured_G(t)
            trans_nablaU * x_G
        end
  end
  
  Ured_G_star
  
  
  
% 2.2 平衡点での Wred_G_star を確認... delta(t,i) -> delta_star(i), E(t,i) -> E_star(i)
  
  Ured_G_inputstar = zeros(sol_size,1);
  trans_nablaU_inputstar = zeros(1,6);
  x_G_inputstar = zeros(6,1);
  Wred_G_inputstar = zeros(sol_size,1);
  
  for t = 1:sol_size
      
        for i = 1:3
              %電気サブシステムGの出力: y_G = -Ei * Σ(j=1,2,3) Ej * Bredij * sinδij
              y_G = 0;
              
              % i 番目の Ured_G を計算する.　　Ured_Gi_star は平衡点での値
              Ured_Gi = (Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i)))) /2;
              
              for j = 1:3
                          
                  Ured_Gi = Ured_Gi + (E_star(i)*E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j))) /2;
                  y_G = y_G - E_star(i) * E_star(j)*Bred(i,j)*sin(delta_star(i)-delta_star(j));
              
              end
              
              % Ured_G = Ured_G1 + Ured_G2 + UredG3
              Ured_G_inputstar(t) = Ured_G_inputstar(t) + Ured_Gi;
              
              % trans_nablaU = [y_G1,y_G2,yG3,Vfield_star(1)/(Xd(1)-Xq(1)),Vfield_star(i)/(Xd(2)-Xq(2)),Vfield_star(3)/(Xd(3)-Xq(3))]
              trans_nablaU_inputstar(i) = y_G;
              trans_nablaU_inputstar(i+3) = Vfield_star(i)/(Xd(i)-Xq(i));
              
              % x_G = [delta(t,1)-delta_star(1); delta(t,2)-delta_star(2); delta(t,3)-delta_star(3); E(t,1)-E_star(1); E(t,2)-E_star(2); E(t,3)-E_star(3)];
              x_G_inputstar(i) = delta_star(i)-delta_star(i);
              x_G_inputstar(i+3) = E_star(i)-E_star(i);

        end
        Wred_G_inputstar(t) = Ured_G_inputstar(t) - Ured_G_star - trans_nablaU_inputstar * x_G_inputstar;

  end

  
  figure;
  plot(t_sol, Ured_G)
  title("U^{red}_G")
  
  figure;
  plot(t_sol, W_F)
  title("W_F")
%{  
  figure;
  plot(t_sol, W_F_star)
  title("W_F *")
%}  
  figure;
  plot(t_sol, Wred_G)
  title("W^{red}_G")
%{  
  figure;
  plot(t_sol, Wred_G_inputstar)
  title("W^{red}_G *")
%}
  figure;
  plot(t_sol,W_F+Wred_G)
  title("W_F + W^{red}_G")

  
  
  % flag_accum_diff == 1 なら [diff( W_F )], [diff( Wred_G )], [diff( W_F + W^{red}_G )] を表示
  if flag_accum_diff == 1
 %{
    dff_F = diff(W_F);
    [sz,~] = size(dff_F);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_F)
    title("diff( W_F )")
    
    dff_Gred = diff(Wred_G);
    [sz,~] = size(dff_Gred);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_Gred)
    title("diff( Wred_G )")

 %}
      

    dff_W_FGred = diff(W_F+Wred_G);
    [sz,~] = size(dff_W_FGred);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t,dff_W_FGred) 
    title("diff( W_F + W^{red}_G )")
    yline(0)
    
    max_diff_W_FGred = max(dff_W_FGred)
    
    
    
  %{
    dff_FGred = diff(W_F+Wred_G);

    [sz,~] = size(dff_FGred);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t, dff_FGred)
    title("diff( W_F + W^{red}_G )")
  %}  
  end



end