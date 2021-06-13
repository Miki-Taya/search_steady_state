%入力を自分で指定し、yを出力とする、Gのみの常微分方程式


global Xd Xq Bred u


%時間区間とYの初期値(Y=[delta(1);delta(2);delta(3);E(1);E(2);E(3)])
tspan = [0 100];    %temp
Y0 = [1;0.1;10;5;3;1];   %temp


%時間tSolに対する数値的に求めた常微分方程式の解YSol
[tSol YSol] = ode45(@generator,tspan,Y0);


%YSolの1~3項は内部電圧E、4~6項は回転子偏角δ
E = [YSol(:,1) YSol(:,2) YSol(:,3)];
delta = [YSol(:,4) YSol(:,5) YSol(:,6)];


%出力yとポテンシャルエネルギー関数U
for t = 1:size(tSol)
  U(t,1) = 0;
  for i = 1:3
    y_temp(t,i) = 0;
    U_temp(t,i) = 0;
    for j = 1:3
      y_temp(t,i) = y_temp(t,i) + E(t,j)*Bred(i,j)*sin(delta(t,i)-delta(t,j));
      U_temp(t,i) = U_temp(t,i) + E(t,j)*Bred(i,j)*cos(delta(t,i)-delta(t,j));
    end
    U(t,1) = U(t,1) + Xd(i)*E(t,i)^2/(Xq(i)*(Xd(i)-Xq(i))) + E(t,i)*U_temp(t,i) ;
  end
  U(t,1) = U(t,1)/2;
end
y = - E.*y_temp;


%グラフ出力：数値的に求めた常微分方程式における時間tSolに対する出力y=[y(1) y(2) y(3)],内部電圧E=[E(1) E(2) E(3)],回転子偏角delta=[delta(1) delta(2) delta(3)]
f1 = figure;
subplot(2,1,1)
plot(tSol,[y delta E])
legend('y1','y2','y3','E1','E2','E3','delta1','delta2','delta3','Location','southwest')
subplot(2,1,2)
plot(tSol,U)
legend('U')

f2 = figure;
subplot(4,1,1)
plot(tSol,y)
ylabel('y')
legend('y1','y2','y3')

subplot(4,1,2)
plot(tSol,delta)
ylabel('delta')
legend('delta1','delta2','delta3')

subplot(4,1,3)
plot(tSol,E)
ylabel('E')
legend('E1','E2','E3')

subplot(4,1,4)
plot(tSol,U)
ylabel('U')
legend('U')


%関数定義：独立変数tと従属変数Y=[delta(1);delta(2);delta(3);E(1);E(2);E(3)]
function dYdt = generator(t,Y)
global Xd Xq Bred u

  %パラメータ設定
  Xd = [1.6;1.4;1.2];   %defined
  Xq = [0.25;0.15;0.15];   %defined
  B = [1 -3 2; 5 -4 -1; 2 1 -2];  %temp.アドミタンス行列Yの虚部であるサセプタンス行列
  Bred = - inv(diag(Xq) - diag(Xq)*B*diag(Xq));
  taud = [5;6;8];   %temp
  Vfield= [1.05;1.1;1];   %temp
  u = [1;4;9];    %temp


  %tの従属変数：delta,E
  delta = [Y(1);Y(2);Y(3)];
  E = [Y(4);Y(5);Y(6)];


  %定義：ddelta/dt
  ddeltadt = u;


  %定義：dE/dt
  for i = 1:3
    for j = 1:3
      A = 0;
      A = A + E(j)*Bred(i,j)*cos(delta(i)-delta(j));
    end
    dEdt(i) = - (Xd(i)/Xq(i)*E(i) - (Xd(i)-Xq(i))*A + Vfield(i)) / taud(i);
  end


  %定義：dY/dt = [ ddelta/dt, dE/dt ]
  dYdt = [ddeltadt;dEdt(1);dEdt(2);dEdt(3)];

end
