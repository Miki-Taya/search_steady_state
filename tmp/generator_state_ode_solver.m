clear

%発電機の内部状態 (δ,E,Δω) の時間応答を表示する
%how to use
%1. ある状態 (δ*,E*,Δω*) から [Pmech*,Vfield*]を計算
%2. [Pmech*,Vfield*] をodeに代入
%3. (δ(0),E(0),Δω(0)) を (δ*,E*,Δω*) からずらしてstart
%4. (δ,E,Δω) の時間応答は (δ*,E*,Δω*) に漸近する


%パラメータ設定
y12 = imag(inv(0.085i));  %1-2間送電線のインピーダンス：z12=0.085j
y23 = imag(inv(0.092i));  %2-3間送電線のインピーダンス：z32=0.092j
Xd = [1.6;1.4;1.2];
Xq = [0.25;0.15;0.15];
B = [y12 -y12 0; -y12 y12+y23 -y23; 0 -y23 y23];  %B：アドミタンス行列Yの虚部であるサセプタンス行列
Bred = - inv(diag(Xq) - diag(Xq)*B*diag(Xq));


%時間区間とgenerator_stateの初期値([delta; E; deltaomega] = [generator_state(1:3); generator_state(4:6); generator_state(7:9)])
tspan = [0 100];
initial_generator_state = [-0.0626;0.5376;2.3497;5.1769;1.4261;2.4998;0;0;0];


%ある内部状態の定常値 (δ*,E*,Δω*)
steady_generator_state = [-0.0626;0.5376;2.3497;5.1769;1.4261;2.4998;0;0;0];
delta_star = steady_generator_state(1:3);
E_star = steady_generator_state(4:6);


%状態 (δ*,E*,Δω*) から [Pmech*,Vfield*]を計算
%[Pmech_star,Vfield_star]を返り値とする関数get_steady_Pmech_Vfieldに、引数(delta_star,E_star,Bred,Xq,Xd)を渡して計算
[Pmech_star,Vfield_star] = get_steady_Pmech_Vfield(delta_star,E_star,Bred,Xq,Xd);


%(t,generator_state,Xd,Xq,Bred,B,Pmech_star,Vfield_star)を引数とする関数generatorをwrapし、引数を(t,generator_state)だけにする
get_dx_nonlinear_ode_wrap = @(t, generator_state) get_dx_nonlinear_ode(t, generator_state, Xd, Xq, Bred, B, Pmech_star, Vfield_star);
%get_dx_nonlinear_ode_wrap(0, steady_generator_state);  %absolutely 0


%時間tSolに対する数値的に求めた常微分方程式の解generator_state_sol
[t_sol generator_state_sol] = ode45(get_dx_nonlinear_ode_wrap, tspan, initial_generator_state);


%generator_state_solの1~3項は回転子偏角δ、4~6項は内部電圧E、7~9項は周波数偏差Δω
delta = generator_state_sol(:,1:3);
E = generator_state_sol(:,4:6);
deltaomega = generator_state_sol(:,7:9);


%ポテンシャルエネルギー関数U[size(t_sol)×1]
for t = 1:size(t_sol)
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


%グラフ出力：数値的に求めた常微分方程式における時間t_solに対する、回転子偏角δ：delta,内部電圧E,周波数偏差Δω：deltaomega,ポテンシャル関数U
f1 = figure;
subplot(4,1,1)
plot(t_sol, delta)
ylabel('delta')
legend('delta1','delta2','delta3')

subplot(4,1,2)
plot(t_sol, E)
ylabel('E')
legend('E1','E2','E3')

subplot(4,1,3)
plot(t_sol, deltaomega)
ylabel('deltaomega')
legend('deltaomega1','deltaomega2','deltaomega3')

subplot(4,1,4)
plot(t_sol, U)
ylabel('U')
legend('U')

%{
f2 = figure;
subplot(2,1,1)
plot(t_sol, [delta E deltaomega])
legend('delta1','delta2','delta3','E1','E2','E3','deltaomega1','deltaomega2','deltaomega3','Location','northwest')
subplot(2,1,2)
plot(t_sol, U)
legend('U')
%}
