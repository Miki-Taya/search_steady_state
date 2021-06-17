function psi = get_psi(steady_generator_state)

    %安定性を定常値を入力

    delta_star = steady_generator_state(1:3);
    E_star = steady_generator_state(7:9);

    %パラメータ設定
    taud = diag([5 6 8]);
    D = diag([2 1.8 2]);
    M = diag([18 13 12]);
    y12 = imag(inv(0.085i));  %1-2間送電線のインピーダンス：z12=0.085j
    y23 = imag(inv(0.092i));  %2-3間送電線のインピーダンス：z32=0.092j
    Xd = [1.6;1.4;1.2];
    Xq = [0.25;0.15;0.15];
    BB = [y12 -y12 0; -y12 y12+y23 -y23; 0 -y23 y23];  %B：アドミタンス行列Yの虚部であるサセプタンス行列
    Bred = - inv(diag(Xq) - diag(Xq)*BB*diag(Xq));
    omega0 = 376.9911;  
    k = zeros(3);
    h = zeros(3);
    L = zeros(3);
    A = zeros(3);
    B = zeros(3);
    C = zeros(3);
    
    
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

    eh = zeros(1,3);
    for i = 1:3
        eh(i) = 2*E_star(i)*h(i,i);
    end

    Xdq = Xd - Xq;
    A = diag(Xdq) * A;
    B = diag(Xdq) * B;
    C = diag(eh) + C;
    
    psi = [zeros(3) omega0*eye(3) zeros(3); -M\L -M\D -M\C; taud\B zeros(3) taud\A];
end
