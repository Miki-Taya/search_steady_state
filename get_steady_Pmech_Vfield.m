
%発電機内部状態[δ,E]から[Pmech_star,Vfield_star]の逆算


function [Pmech_star, Vfield_star] = get_steady_Pmech_Vfield(delta, E, Bred, Xq, Xd)

  for i = 1:3
    sigma_sin = 0;  %f/(-E)
    sigma_cos = 0;  %-g
    for j = 1:3
      sigma_sin = sigma_sin + E(j)*Bred(i,j)*sin(delta(i)-delta(j));
      sigma_cos = sigma_cos + E(j)*Bred(i,j)*cos(delta(i)-delta(j));
    end
    Pmech_star(i) = - E(i)*sigma_sin;
    Vfield_star(i) = Xd(i)/Xq(i)*E(i) + (Xd(i)-Xq(i))*sigma_cos;
  end

end
