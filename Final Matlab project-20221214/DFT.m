function U=DFT(u,M)
    Nx=length(u);
    omega=2*pi/M*(-M/2+1:M/2);
    U=zeros(M,1);
    for n=1:M
        omega_n=omega(n);
        for k_x=1:Nx
            U(n)=U(n)+u(k_x)*exp(-i*omega_n*k_x);
        end
    end
end