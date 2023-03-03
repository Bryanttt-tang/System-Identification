function plot_mag(u)
    Ts=0.1; M3=2047;
    u1=reshape(u,M3,[]);
    u_final=mean(u_1,2);
    U=fft(u_final);
    omega=2*pi/(Ts*M3)*(36:175)
    plot(omega,20*log10(abs(U(36:175)))
end