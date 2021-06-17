function judge_steady(generator_state)
%安定性を定常値を入力

%安定
%delta_star = [1.1;1.6;1.4];
%E_star = [1;2.5;2.1];

%安定
%steady_delta = [0.15708;0;-0.15708];
%steady_E = [3.2225;3.2225;3.2225];

%不安定
%steady_delta = [0.942478;0;-1.41372];
%steady_E = [3.2225;3.2225;3.2225];

delta = generator_state(1:3);
E = generator_state(7:9);

%パラメータ設定
taud = diag([5.1400, 5.9000, 8.9700]);
D = diag([2, 2, 2]);
M = diag([100, 18, 12]);
Xq = [0.9360;0.9110;0.6670];
Xd = [1.5690;1.6510;1.2200];
BB = [-6.1331,1.4914,1.6779; 1.4914,-5.9131,2.2693; 1.6779,2.2693,-5.6149];  %BB：アドミタンス行列Yの虚部であるサセプタンス行列
Bred = - inv(diag(Xq) - diag(Xq)*BB*diag(Xq));
omega0 = 376.9911;  

for i = 1:3
    for j = 1:3
        k(i,j) = -Bred(i,j)*cos(delta(i) - delta(j));
        h(i,j) = -Bred(i,j)*sin(delta(i) - delta(j));
    end
end

 

for i = 1:3
    for j = 1:3
        if i == j
            Ek = 0; Eh = 0;
            for q = 1:3
                if i == q
                    continue
                end
                
                Ek = Ek + E(q) * k(i,q);
                Eh = Eh + E(q) * h(i,q);                
            end
            
            L(i,j) = E(i) * Ek; 
            A(i,j) = k(i,i) - Xd(i)/(Xq(i)*(Xd(i) - Xq(i)));
            B(i,j) = - Eh;
            C(i,j) = Eh; 
            
        else
            
             L(i,j) = - E(i) * E(j) * k(i,j); 
             A(i,j) = k(i,j);
             B(i,j) = E(j) * h(i,j); 
             C(i,j) = E(i) * h(i,j);
             
        end
    end
end

for i = 1:3
    eh(i) = 2*E(i)*h(i,i);
end

Xdq = Xd - Xq;
A = diag(Xdq) * A;
B = diag(Xdq) * B;
C = diag(eh) + C;

psi = [zeros(3) omega0*eye(3) zeros(3); -M\L -M\D -M\C; taud\B zeros(3) taud\A];

lamdapsi = eig(psi);
lamdaA = eig(A);
lamdaB = eig(B);
lamdaL = eig(L);

sum(L,2);
sum(B,2);


%必要な条件
% linear:[psi is steady](線形システムの漸近安定性) and [matrix:L * span{1} = 0, matrix:B * span{1} = 0]
% feedback Gの受動性:[A is steady] and [Lo is symmetric positive semi-definite](受動送電条件1,3)
% non-linear: 蓄積関数 W が半正定関数

% psi や A の安定性を判定
if all(real(lamdapsi) < 0)
    disp('Matrix: psi is steady.');
else
    disp('Matrix: psi is not steady.');
end

if all(real(lamdaA) < 0)
    disp('Matrix: A is steady.');
else
    disp('Matrix: A is not steady.');
end

% L や B のカーネルが 1 ( matrix:L * span{1} = 0 )かどうかを判定
if all(abs(sum(B,2)) < 10^(-12)) % e-15 程度の誤差がある
    disp('Kernel of Matrix:B is span{1}.');
else
    disp('Kernel of Matrix:B is not span{1}.');
end

if all(abs(sum(L,2)) < 10^(-12))
    disp('Kernel of Matrix:L is span{1}.');
else
    disp('Kernel of Matrix:L is not span{1}.');
end

%Loが対称正定値行列かどうかを判定（semi-definiteより厳しい条件）
Lo = L - C*(A\B);

try chol(Lo); %対称正定値行列を対称部分と上三角部分だけ使用して表す関数 chol を呼ぶ
    disp('Matrix: Lo=L-C*inv(A)*B　is symmetric positive definite.') % disp は値を表示 
catch ME % try error である MException の略
    if all(eig(Lo) >= 0) & (Lo - transpose(Lo) < 10^(-10))
      disp('Matrix: Lo=L-C*inv(A)*B　is symmetric positive semi-definite')
    else
      disp('Matrix: Lo=L-C*inv(A)*B　is not symmetric positive semi-definite')
    end
end





