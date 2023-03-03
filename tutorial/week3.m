N=4096; Ts=1; %sampling time
%%
e_k=randn(N,1);
E_n=fourier_trans(e_k,length(e_k));
% periodog = periodogram(e_k) ;
omega = [] ;% non-negative frequency
for i = 1:N
    omega(i,1) = 2*pi*(i-1)/N ;
end
idx = find(omega > 0 & omega < pi); % positive frequencies
peri=abs(E_n(idx)).^2/N;
figure(1)
loglog(omega(idx),peri);
%%
z = tf('z',1);
P_z = (z + 0.5) / ((z+0.5)*(z-0.5)^2+ 0.5);
t_s=[];
for i = 1:N
    t_s(i,1) = Ts*(i-1);
end
w=lsim(P_z,e_k,t_s);
W=fft(w);
peri_w=abs(W(idx)).^2/N;
%%
[mag,phase,wout] = bode(P_z,omega);
mag=squeeze(mag(1,1,:));
error= peri_w-abs(mag(idx));
figure(2)
loglog(omega(idx),peri_w);
hold on
loglog(omega(idx),abs(mag(idx)));
hold on
loglog(omega(idx),error);
legend('peri','bode','error')
%%
function U_n=fourier_trans(u,N)
U_n=fft(u,N);
end