%安定性を定常値を入力
steady_delta = [0.1;1;-2];
steady_E = [-3;1;2];

%パラメータ設定
taud = diag([5;6;8]);
D = diag([2 1.8 2]);
M = diag([18 13 12]);
y12 = imag(inv(0.085i));  %1-2間送電線のインピーダンス：z12=0.085j
y23 = imag(inv(0.092i));  %2-3間送電線のインピーダンス：z32=0.092j
Xd = [1.6;1.4;1.2];
Xq = [0.25;0.15;0.15];
B = [y12 -y12 0; -y12 y12+y23 -y23; 0 -y23 y23];  %B：アドミタンス行列Yの虚部であるサセプタンス行列
Bred = - inv(diag(Xq) - diag(Xq)*B*diag(Xq));
omega0 = 376.9911;  

for i = 1:3
    for j = 1:3
        k(i,j) = -Bred(i,j)*cos(steady_delta(i) - steady_delta(j));
        h(i,j) = -Bred(i,j)*sin(steady_delta(i) - steady_delta(j));
    end
end

 

for i = 1:3
    for j = 1:3
        if i == j
            Ek = 0; Eh = 0;
            for j = 1:3
                Ek = Ek + steady_E(j) * k(i,j);
                Eh = Eh + steady_E(j) * h(i,j);                
            end
            
            L(i,j) = steady_E(i) * Ek; 
            A(i,j) = k(i,i) - Xd(i)/(Xq(i)*(Xd(i) - Xq(i)));
            B(i,j) = - Eh; 
            C(i,j) = Eh; 
            
        else
            
             L(i,j) = steady_E(i) * steady_E(j) * k(i,j); 
             A(i,j) = k(i,i);
             B(i,j) = steady_E(j) * h(i,j); 
             C(i,j) = steady_E(i) * h(i,j);
             
        end
    end
end



psi = [zeros(3) omega0*eye(3) zeros(3); -inv(M)*L -inv(M)*D -inv(M)*C; inv(taud)*B zeros(3) inv(taud)*A];

lamda = eig(psi)