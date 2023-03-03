
clc; clear all; close all;
[p3_Tur, p3_ID, p3_beta, p3_alpha, p3_Rhat, p3_Ghat,p3_r,p3_y,p3_wl,p3_wu,p3_a,R] = HS2022_SysID_final_p3_21957907();

% G_s = tf(conv([1 2*zeta_z*omega_z omega_z^2], [5000]), conv([1 2*zeta_p*omega_p omega_p^2], conv([1 50], [1 200])));
% G_dz = c2d(G_s, T_s, 'zoh');
Ts=1;
% C_dz = tf([1.25 -0.75], [1 -1], 1);
% C=tf([1.25 -0.75], [1 -1]);
N=size(p3_r,1);
% figure(88)
% r1=p3_r(:,1); r2=p3_r(:,2); r3=p3_r(:,3);
% omega_k= 2*pi/(Ts*N)*(1:(N/2));
% range=find(omega_k>=p3_wl & omega_k<=p3_wu);
% range_plus=[range,max(range)+1]; % include wu
% plot(range_plus,r1(range_plus)); hold on;
% plot(range_plus,r2(range_plus)); hold on;
% plot(range_plus,r3(range_plus)); hold on;
% legend('1','2','3') 

%%
  a=p3_a(p3_ID); r=p3_r(:,p3_ID); y=p3_y(:,p3_ID); k=N/2;
  omega_full= 2*pi/(Ts*N)*(1:N);
     C_0 = tf(a, [1 -1], Ts);
     u=r-lsim(C_0,y);
     U1=fft(u); R1=fft(r); Y1=fft(y);
     est_1=U1(2:end)./R1(2:end);
     Tyr=Y1(2:end)./R1(2:end);
     
     C_fre=squeeze(freqresp(C_0,omega_full));
     est_2=1-C_fre(1:end-1).*Tyr;
