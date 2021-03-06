function plot_FG_accum_func(steady_generator_state, delta, deltaomega, E, B_sus ,Bred,Y, Xd, Xq, Vfield_star, omega0, M, t_sol, flag_accum_diff)
  
  delta_star = steady_generator_state(1:3);
  deltaomega_star = transpose(steady_generator_state(4:6));
  E_star = steady_generator_state(7:9);
  
  [sol_size,~] = size(t_sol);

%---------------------------------------------------------------------------
  
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
  
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
  %{ 
%G：非線形微分代数方程式から発電機バスのクロン縮約をして、非線形常微分方程式系で考えたポテンシャルエネルギー関数(Ured_G)と蓄積関数(Wred_G)

% 1.1 matrix を for文で表した　こっちはよくない


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
  
  %}
  

  %{  
% 1.2 平衡点を代入したら W は 0 にならなければならない。そのテスト。for ver.


  
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
         
              % Ured_G に項を加えていく.　　Ured_Gi_star は平衡点での値
              Ured_G(t) = Ured_G(t) + (Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i)))) /2;
              
              if t == 1
                  
                 % Ured_Gi_star に項を加えていく。平衡点での値(star)は t によらないため、t=1 のときだけ 
                 Ured_G_star = Ured_G_star + (Xd(i)*E_star(i)^2/(Xq(i)*(Xd(i)-Xq(i)))) /2;
                 
                 % trans_nablaU(4~6) = Vfield_star(i)/(Xd(i)-Xq(i)) を代入
                 trans_nablaU(i+3) = Vfield_star(i)/(Xd(i)-Xq(i));
                 
              end
                 
              for j = 1:3
                  
                  % 平衡点での値(star)は t によらないため、t == 1 のときだけ計算
                  if t == 1
                      
                        Ured_G_star = Ured_G_star + (E_star(i)*E_star(j)*Bred(i,j)*cos(delta_star(i)-delta_star(j))) /2;

                        %電気サブシステムGの出力: y_G = -Ei * Σ(j=1,2,3) Ej * Bredij * sinδij                                                  
                        % trans_nablaU(1~3) に y_G_star(i) を代入
                        trans_nablaU(i) = trans_nablaU(i) - E_star(i)*E_star(j)*Bred(i,j)*sin(delta_star(i)-delta_star(j));                            

                  end
                       
                  Ured_G(t) = Ured_G(t) + (E(t,i)*E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j))) /2;
              
              end
               
              % x_G = [delta(t,1)-delta_star(1); delta(t,2)-delta_star(2); delta(t,3)-delta_star(3); E(t,1)-E_star(1); E(t,2)-E_star(2); E(t,3)-E_star(3)];
              x_G(i) = delta(t,i)-delta_star(i);
              x_G(i+3) = E(t,i)-E_star(i);

        end
        Wred_G(t) = Ured_G(t) - Ured_G_star - trans_nablaU * x_G;
        %{
        %Eの初期誤差0のときに応答が振動してしまうｋら、Ured_G , trans_nablaU * x_G の値を調べる
        if mod(t,300) == 0
            Ured_G(t)
            trans_nablaU * x_G
        end
        %}
  end
    
  
%{  
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
%}
  
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
  
% 微分代数方程式系の U_G を計算する
  
% まずは 電圧フェーザ:V (3,sol_size)(V(t) は縦ベクトル) から V_abs, V_arg を求める

E = transpose(E);  % 行列計算のために E(t)を縦ベクトルにする (3*sol_size )
delta = transpose(delta);  % 行列計算のために delta(t)を縦ベクトルにする(3*sol_size )
V = zeros(3,sol_size); % (3*sol_size )


for t = 1:sol_size

    V(:,t) = (diag(1./(1j*Xq)) + Y) \ diag(exp(1j*delta(:,t))./(1j*Xq)) * E(:,t);
    
end

V_star = (diag(1./(1j*Xq)) + Y) \ (diag(exp(1j*delta_star)./(1j*Xq))) * E_star;

V_abs = abs(V);
V_arg = angle(V);
V_abs_star = abs(V_star);
V_arg_star = angle(V_star);

%{
% V_star を使って、微分代数方程式の右辺が 0 になるかを見る

tauE_star = zeros(3,1);

for i = 1:3
    tauE_star(i) = -Xd(i)/Xq(i)*E_star(i) + (Xd(i)/Xq(i)-1)*V_abs_star(i)...
                   *cos(delta_star(i)-V_arg_star(i)) + Vfield_star(i);
end

%}
% V を使って U_G を求める

U_G = zeros(sol_size,1);
U_G_star = 0;

for t = 1:sol_size
    
    for i = 1:3
        
        if t == 1
            U_G_star = U_G_star + Xd(i)*E_star(i)^2/(2*Xq(i)*(Xd(i)-Xq(i)))...
                - E_star(i)*V_abs_star(i)/Xq(i)*cos(delta_star(i)-V_arg_star(i))...
                + V_abs_star(i)^2/(2*Xq(i)) - B_sus(i,i)*V_abs_star(i)^2/2; 
        end
        
        U_G(t) = U_G(t) + Xd(i)*E(i,t)^2/(2*Xq(i)*(Xd(i)-Xq(i)))...
            - E(i,t)*V_abs(i,t)/Xq(i)*cos(delta(i,t)-V_arg(i,t))...
            + V_abs(i,t)^2/(2*Xq(i)) - B_sus(i,i)*V_abs(i,t)^2/2;
        
        for j = 1:3    
            
            if j ~= i      
                
                if t == 1
                    U_G_star = U_G_star - 1/2 * B_sus(i,j)*V_abs_star(i)*V_abs_star(j)*cos(V_arg_star(i)-V_arg_star(j));
                end
                
                U_G(t) = U_G(t) - 1/2 * B_sus(i,j)*V_abs(i,t)*V_abs(j,t)*cos(V_arg(i,t)-V_arg(j,t));
            
            end       
        end        
    end   
end
E = transpose(E);  % E(t)を横ベクトルに戻す  (sol_size*3)
delta = transpose(delta);  % delta(t)を横ベクトルに戻す  (sol_size*3)
  
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% 微分代数方程式系
% U_G から W_G を求める

  trans_nablaU = zeros(1,6);
  x_G = zeros(6,1);
  W_G = zeros(sol_size,1);


  for t = 1:sol_size
      
        for i = 1:3
            
                if t == 1
                      %電気サブシステムGの出力: y_G_star を trans_nablaU(1~3)に代入
                      % y_G = E(i)*V_abs(i)/Xq(i)*sin(delta(i) - V_arg(i))
                      % trans_nablaU = [y_G1_star,y_G2_star,yG3_star,Vfield_star(1)/(Xd(1)-Xq(1)),Vfield_star(i)/(Xd(2)-Xq(2)),Vfield_star(3)/(Xd(3)-Xq(3))]
                      trans_nablaU(i) = E_star(i)*V_abs_star(i)/Xq(i)*sin(delta_star(i)-V_arg_star(i));              
                      trans_nablaU(i+3) = Vfield_star(i)/(Xd(i)-Xq(i));
                end
                
                % x_G = [delta(t,1)-delta_star(1); delta(t,2)-delta_star(2); delta(t,3)-delta_star(3); E(t,1)-E_star(1); E(t,2)-E_star(2); E(t,3)-E_star(3)];
                x_G(i) = delta(t,i)-delta_star(i);
                x_G(i+3) = E(t,i)-E_star(i);

        end

        W_G(t) = U_G(t) - U_G_star - trans_nablaU * x_G;

  end
%---------------------------------------------------------------------------
%--------------------------------------------------------------------------- 
% DAE と ODE の差分

  figure;
  plot(t_sol,U_G-Ured_G)
  title("U_G - Ured_G")
  
  figure;
  plot(t_sol,W_G - Wred_G)
  title("W_G - Wred_G") 
  
% W_F
  figure;
  plot(t_sol, W_F)
  title("W_F")
 
  
% 同じグラフ内で表示して比較
  figure;
  subplot(1,3,1)
  plot(t_sol, [Ured_G, U_G])
  title("U^{red}_G,U_G")
  legend("U^{red}_G","U_G")
  
  subplot(1,3,2)
  plot(t_sol, [W_F, Wred_G, W_G])
  title("W_F, W^{red}_G, W_G")  
  legend("W_F","W^{red}_G","W_G")

  subplot(1,3,3)
  plot(t_sol,[W_F+Wred_G, W_F+W_G])
  title("W_F + W^{red}_G, W_F+W_G")  
  legend("W_F + W^{red}_G","W_F+W_G")


  %別々に表示
  figure;
  subplot(2,3,1)
  plot(t_sol, Ured_G)
  title("U^{red}_G")
  
  subplot(2,3,2)
  plot(t_sol, Wred_G)
  title("W^{red}_G")
  
  subplot(2,3,3)
  plot(t_sol,W_F+Wred_G)
  title("W_F + W^{red}_G")
  
  subplot(2,3,4)
  plot(t_sol, U_G)
  title("U_G")
  
  subplot(2,3,5)
  plot(t_sol, W_G)
  title("W_G")
  
  subplot(2,3,6)
  plot(t_sol,W_F+W_G)
  title("W_F + W_G") 
  

  %別々に表示 W_F,W_G,W_F+W_G;W_F,Wred_G,W_F+Wred_G
  figure;
  subplot(2,3,1)
  plot(t_sol, W_F)
  
  subplot(2,3,2)
  plot(t_sol, W_G)
  
  subplot(2,3,3)
  plot(t_sol,W_F+W_G, "LineWidth", 2)
  
  subplot(2,3,4)
  plot(t_sol, W_F)
  
  subplot(2,3,5)
  plot(t_sol, Wred_G)
  
  subplot(2,3,6)
  plot(t_sol,W_F+Wred_G, "LineWidth", 2)

  
%{  
  figure;
  plot(t_sol, W_F_star)
  title("W_F *")
%}  

  

%{  
  figure;
  plot(t_sol, Wred_G_inputstar)
  title("W^{red}_G *")
%}

  


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------  
  
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
      
    %一緒に表示
    dff_W_FGred = diff(W_F+Wred_G);
    dff_W_FG = diff(W_F+W_G);
    [sz,~] = size(dff_W_FGred);
    t = transpose(linspace(0,100,sz));
    figure;
    plot(t,[dff_W_FGred, dff_W_FG]) 
    title("diff( W_F + W^{red}_G ), diff( W_F + W_G )")
    yline(0)    
    legend("diff( W_F + W^{red}_G )","diff( W_F + W_G )")
    

    %別々に表示
    dff_W_FGred = diff(W_F+Wred_G);
    [sz,~] = size(dff_W_FGred);
    t = transpose(linspace(0,100,sz));
    figure;
    subplot(1,2,1)
    plot(t,dff_W_FGred) 
    title("diff( W_F + W^{red}_G )")
    yline(0)
    
    dff_W_FG = diff(W_F+W_G);
    [sz,~] = size(dff_W_FG);
    t = transpose(linspace(0,100,sz));
    subplot(1,2,2)
    plot(t,dff_W_FG) 
    title("diff( W_F + W_G )")
    yline(0)

      
% dff_W_FGred は常に負であってほしい。max を取ってそれを調べる
    max_diff_W_FGred = max(dff_W_FGred)
    max_diff_W_FG = max(dff_W_FG)

  end

end