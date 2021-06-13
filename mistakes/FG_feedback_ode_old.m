clear
global Xd Xq Bred B


%時間区間とYの初期値(Y=[delta(1);delta(2);delta(3);E(1);E(2);E(3);deltaomega(1);deltaomega(2);deltaomega(3))]
tspan = [0 100];    %temp
Y0 = [-0.0626;0.5376;2.3497;5.1769;1.4261;2.4998;0;0;0];   %定常値：Y0_star = [-0.0626;0.5376;2.3497;5.1769;1.4261;2.4998;0;0;0]


%時間tSolに対する数値的に求めた常微分方程式の解YSol
[tSol YSol] = ode45(@generator,tspan,Y0);


%YSolの1~3項は回転子偏角δ、4~6項は内部電圧E、7~9項は周波数偏差Δω
delta = [YSol(:,1) YSol(:,2) YSol(:,3)];
E = [YSol(:,4) YSol(:,5) YSol(:,6)];
deltaomega = [YSol(:,7) YSol(:,8) YSol(:,9)];


%ポテンシャルエネルギー関数U[tSol×1]
for t = 1:size(tSol)
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


%グラフ出力：数値的に求めた常微分方程式における時間tSolに対する、回転子偏角δ：delta=[delta(1) delta(2) delta(3)],内部電圧E=[E(1) E(2) E(3)],周波数偏差Δω：deltaomega=[deltaomega(1);deltaomega(2);deltaomega(3)]
f1 = figure;
subplot(4,1,1)
plot(tSol,delta)
ylabel('delta')
legend('delta1','delta2','delta3')

subplot(4,1,2)
plot(tSol,E)
ylabel('E')
legend('E1','E2','E3')

subplot(4,1,3)
plot(tSol,deltaomega)
ylabel('deltaomega')
legend('deltaomega1','deltaomega2','deltaomega3')

subplot(4,1,4)
plot(tSol,U)
ylabel('U')
legend('U')

%{
f2 = figure;
subplot(2,1,1)
plot(tSol,[delta E deltaomega])
legend('delta1','delta2','delta3','E1','E2','E3','deltaomega1','deltaomega2','deltaomega3','Location','northwest')
subplot(2,1,2)
plot(tSol,U)
legend('U')
%}

%関数定義：独立変数tと従属変数Y=[delta(1);delta(2);delta(3);E(1);E(2);E(3);deltaomega(1);deltaomega(2);deltaomega(3)]
function dYdt = generator(t,Y)
global Xd Xq Bred B

  %パラメータ設定
  Xd = [1.6;1.4;1.2];
  Xq = [0.25;0.15;0.15];
  y12 = imag(inv(0.085i));  %1-2間送電線のインピーダンス：z12=0.085j
  y23 = imag(inv(0.092i));  %2-3間送電線のインピーダンス：z32=0.092j
  B = [y12 -y12 0; -y12 y12+y23 -y23; 0 -y23 y23];  %B：アドミタンス行列Yの虚部であるサセプタンス行列
  Bred = - inv(diag(Xq) - diag(Xq)*B*diag(Xq));
  taud = [5;6;8];
  %Vfield_star = [32.6983;-1.1093;29.3787];
  Vfield_star = [26.0668;12.8409;13.0063];   %定常値：Vfield_star = [32.6983;-1.1093;29.3787]or[26.0668;12.8409;13.0063]
  Pmech_star = [-1.6967 0.8984 0.7984];
  D = [2 1.8 2];
  M = [18 13 12];
  omega0 = 376.9911;   %系統周波数：60[Htz]*2pi


  %tの従属変数「δ：delta, E, Δω：deltaomega」
  delta = [Y(1);Y(2);Y(3)];
  E = [Y(4);Y(5);Y(6)];
  deltaomega = [Y(7);Y(8);Y(9)];


  %定義：ddelta/dt,dE/dt,ddeltaomega/dt
  for i = 1:3
    temp_E = 0;
    temp_deltaomega = 0;
    for j = 1:3
      temp_E = temp_E + E(j)*Bred(i,j)*cos(delta(i)-delta(j));
      temp_deltaomega = temp_deltaomega + E(j)*Bred(i,j)*sin(delta(i)-delta(j));
    end
    ddeltadt(i) = omega0*deltaomega(i);
    dEdt(i) = (-Xd(i)/Xq(i)*E(i) - (Xd(i)-Xq(i))*temp_E + Vfield_star(i)) / taud(i);
    ddeltaomegadt(i) = (E(i)*temp_deltaomega - D(i)*deltaomega(i) + Pmech_star(i)) / M(i);
  end


  %定義：dY/dt = [ ddelta/dt; dE/dt; ddeltaomega/dt]
  dYdt = [ddeltadt(1);ddeltadt(2);ddeltadt(3);dEdt(1);dEdt(2);dEdt(3);ddeltaomegadt(1);ddeltaomegadt(2);ddeltaomegadt(3)];

end
