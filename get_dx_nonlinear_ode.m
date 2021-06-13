
%KronReductionされた非線形常微分方程式
%独立変数t
%従属変数x: generator_state = [delta; E; deltaomega] (9×1行列)
%mainで代入した定数 (Xd, Xq, Bred, B, Pmech_star, Vfield_star) を引数として渡す
%dx = [ ddelta; dE; ddeltaomega]（9×1行列） を返す


function dx = get_dx_nonlinear_ode(t, x, Xd, Xq, Bred, Pmech_star, Vfield_star)

  %パラメータ設定
  taud = [5.1400, 5.9000, 8.9700];
  D = [2, 2, 2];
  M = [100, 18, 12];
  omega0 = 376.9911;   %系統周波数：60[Htz]*2pi

  %tの従属変数「δ：delta, E, Δω：deltaomega」
  delta = x(1:3);
  deltaomega = x(4:6);
  E = x(7:9);

  ddelta = zeros(1,3);
  ddeltaomega = zeros(1,3);
  dE = zeros(1,3);

  %定義：ddelta/dt, dE/dt, ddeltaomega/dt
  for i = 1:3
    sigma_cos = 0;
    sigma_sin = 0;
    for j = 1:3
      sigma_cos = sigma_cos + E(j)*Bred(i,j)*cos(delta(i)-delta(j));
      sigma_sin = sigma_sin + E(j)*Bred(i,j)*sin(delta(i)-delta(j));
    end
    ddelta(i) = omega0*deltaomega(i);
    ddeltaomega(i) = (E(i)*sigma_sin - D(i)*deltaomega(i) + Pmech_star(i)) / M(i);
    dE(i) = (-Xd(i)/Xq(i)*E(i) - (Xd(i)-Xq(i))*sigma_cos + Vfield_star(i)) / taud(i);
  end


  %定義：dx/dt = [ ddelta/dt; dE/dt; ddeltaomega/dt]
  dx = transpose([ddelta, ddeltaomega, dE]);

end
