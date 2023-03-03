%% Transfer function
Ts=1; z=tf('z',Ts);
G_z=(0.0565*z^-1 + 0.1013*z^-2) / (0.8864*z^-4-1.3161*z^-3+1.5894*z^-2-1.4183*z^-1+1);
%%
N=200;
u1=idinput(N); u2=rand(N,1);
y1= lsim(G_z,u1)+ randn(N,1)* sqrt(0.01); % noisy output
y2= lsim(G_z,u2)+ randn(N,1)* sqrt(0.01);
t_vec=(0:N-1)';
%Autocorrelation assuming periodic input
figure(6)
tau=-N/2+1:N/2;
R_u=zeros(N,1);
for t=1:length(tau)
    for k=1:N
        idx=k-tau(t);
        if idx<1
            idx=k-tau(t)+N;
        elseif idx>N
            idx=k-tau(t)-N;
        end
        R_u(t,1)=R_u(t,1)+1/N*u1(k)*u1(idx);
    end
end
plot(tau,R_u)

% plot time-domain input and output
figure(1)
plot(t_vec,u1); hold on; plot(t_vec,y1);
xlim([0,N]); legend('u','y'); grid on
title('Using random binary input signal')
figure(2)
plot(t_vec,u2); hold on; plot(t_vec,y2);
xlim([0,N]); legend('u','y'); grid on
title('Using random uniform input signal')
%%  b)
Y1=fft(y1); U1=fft(u1); G1_N=Y1(1:N/2+1)./U1(1:N/2+1); %only positive frequencies
Y2=fft(y2); U2=fft(u2); G2_N=Y2(1:N/2+1)./U2(1:N/2+1); %only positive frequencies
omega=2*pi/(Ts*N)*(0:N/2); % only positive omega
% frequency response of true system (Bode or freqresp)
[mag,phase,wout] = bode(G_z,omega);
G_mag=squeeze(mag(1,1,:));
%plot ETFT, true system and magnitue of errors
figure;
subplot(2,1,1)
loglog(omega,abs(G1_N),'b--'); hold on
loglog(omega,abs(G2_N),'r--'); hold on
loglog(omega,abs(G_mag),'green'); hold off
xlim([0.2,max(omega)]); ylim([1e-2,1e2]);
title('ETFE vs true plant')
xlabel('Frequency (rad/s)'); ylabel('Magnitude');
legend('ETFE random binary','ETFE random uniform','True plant' )

subplot(2,1,2)
loglog(omega,abs(G1_N-G_mag),'b--'); hold on
loglog(omega,abs(G2_N-G_mag),'r--'); hold off
xlim([0.2,max(omega)]); ylim([1e-2,1e2]);
title('Magnitude of errors')
xlabel('Frequency (rad/s)'); ylabel('Magnitude');
legend('Error random binary','Error random uniform')
%% c)
% 1) long input
u3=idinput(5*N); y3=lsim(G_z,u3)+randn(5*N,1)*0.1;
Y3=fft(y3); U3=fft(u3); G3_N=Y3(1:5*N/2+1)./U3(1:5*N/2+1); %only positive frequencies
omega2= 2*pi/(Ts*5*N)*(0:5*N/2); %positive frequencies with more frequency points on a unit circle
[mag2,phase2,wout2] = bode(G_z,omega2);
G_mag2=squeeze(mag2(1,1,:)); %true plant frequency response

% 2)repeat the same inputs 5 times. The true response is the same, equivalently
% averaging the noise
u4=u1; y4= lsim(G_z,u4)+mean(randn(N,5)*0.1,2); % average noise along the column
Y4=fft(y4); U4=fft(u4); G4_N=Y4(1:N/2+1)./U4(1:N/2+1);

% 3) concat the same inputs for 5 times to get a periodic signal
u5=u1; % tiling u1 5 times into a long vector
yp=lsim(G_z,repmat(u1,5,1))+randn(5*N,1)*0.1;
%then discard the first period and average the last four periods along the
%column direction
yp_1=reshape(yp(N+1:end),N,[]);
y5=mean(yp_1,2);
Y5=fft(y5); U5=fft(u5); G5_N=Y5(1:N/2+1)./U5(1:N/2+1);
%plot ETFE, true system and errors
figure(5);
subplot(2,1,1)
loglog(omega2,abs(G3_N),'b--'); hold on
loglog(omega,abs(G4_N),'r--'); hold on
loglog(omega,abs(G5_N),'k'); hold on
loglog(omega2,abs(G_mag2),'g'); hold off
xlim([0.2,max(omega)]); ylim([1e-2,1e2]);
title('ETFE vs true plant')
xlabel('Frequency (rad/s)'); ylabel('Magnitude');
legend('ETFE 1','ETFE 2','ETFE 3','True plant' )

subplot(2,1,2)
loglog(omega2,abs(G3_N-G_mag2),'b--'); hold on
loglog(omega,abs(G4_N-G_mag),'r--'); hold on
loglog(omega,abs(G5_N-G_mag),'k'); hold off
xlim([0.2,max(omega)]); ylim([1e-3,1e2]);
title('Magnitude of errors')
xlabel('Frequency (rad/s)'); ylabel('Magnitude');
legend('Error 1','Error 2','Error 3')