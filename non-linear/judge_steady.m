clear

%安定性を定常値を入力

%delta_star = [1.1;1.6;1.4];
%E_star = [1;2.5;2.1];

delta_star = [-1.25664;0;0];
E_star = [3.2225;3.2225;3.2225];

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
        k(i,j) = -Bred(i,j)*cos(delta_star(i) - delta_star(j));
        h(i,j) = -Bred(i,j)*sin(delta_star(i) - delta_star(j));
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
                
                Ek = Ek + E_star(q) * k(i,q);
                Eh = Eh + E_star(q) * h(i,q);                
            end
            
            L(i,j) = E_star(i) * Ek; 
            A(i,j) = k(i,i) - Xd(i)/(Xq(i)*(Xd(i) - Xq(i)));
            B(i,j) = - Eh;
            C(i,j) = Eh; 
            
        else
            
             L(i,j) = - E_star(i) * E_star(j) * k(i,j); 
             A(i,j) = k(i,j);
             B(i,j) = E_star(j) * h(i,j); 
             C(i,j) = E_star(i) * h(i,j);
             
        end
    end
end

for i = 1:3
    eh(i) = 2*E_star(i)*h(i,i);
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



%必要な条件
% linear:[psi is steady](線形システムの漸近安定性) and [L and B have zero eigenvalue]
% feedback:[A is steady] and [Lo is symmetric positive semi-definite](受動送電条件1,3)
% non-linear: 蓄積関数 W が半正定関数

% psi や A の安定性や L や B がゼロ固有値を持つかどうかを判定
if all(real(lamdapsi) < 0)
    disp('Matrix: psi is steady.');
end

if all(real(lamdaA) < 0)
    disp('Matrix: A is steady.');
end

if any(abs(lamdaB) < 10^(-15)) % e-17 程度の誤差がある
    disp('Matrix: B has zero eigenvalue.');
end

if any(abs(lamdaL) < 10^(-15))
    disp('Matrix: L has zero eigenvalue.');
end

%Loが対称正定値行列かどうかを判定（semi-definiteより厳しい条件）
Lo = L - C*inv(A)*B;

try chol(Lo); %対称正定値行列を対称部分と上三角部分だけ使用して表す関数 chol を呼ぶ
    disp('Matrix: Lo=L-C*inv(A)*B　is symmetric positive definite.') % disp は値を表示 
catch ME % try error である MException の略
    disp('Matrix: Lo=L-C*inv(A)*B　is not symmetric positive definite')
end





