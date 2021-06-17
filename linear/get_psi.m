function psi = get_psi(steady_generator_state)

    %安定性を定常値を入力

    delta_star = steady_generator_state(1:3);
    E_star = steady_generator_state(7:9);

    %パラメータ設定
    taud = diag([5.1400, 5.9000, 8.9700]);
    D = diag([2, 2, 2]);
    M = diag([100, 18, 12]);
    Xq = [0.9360;0.9110;0.6670];
    Xd = [1.5690;1.6510;1.2200];
    BB = [-6.1331,1.4914,1.6779; 1.4914,-5.9131,2.2693; 1.6779,2.2693,-5.6149];  %BB：アドミタンス行列Yの虚部であるサセプタンス行列
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
