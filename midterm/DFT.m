function U=DFT(u,M)
    N=length(u);
    omega=2*pi/M*(-(M-1)/2:(M-1)/2);
    U=zeros(M,1)
    for n=1:M
        omega_n=omega(n);
        for lamda=1:N
            U(n)=U(n)+u(lamda)*exp(-j*omega_n*lamda);
        end
    end
end