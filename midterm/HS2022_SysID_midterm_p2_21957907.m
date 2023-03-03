function [p21_ID, p21_omega, p21_G_hat, p22_ID, p22_omega, p22_phi_u, p22_phi_yu, p22_G_hat, p23_ID, p23_omega, p23_G_hat] = HS2022_SysID_midterm_p2_21957907()

%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[p2_data] = HS2022_SysID_midterm_p2_GenerateData(LegiNumber);
%extract each experiment data
u1=p2_data(1).p2_u;u2=p2_data(2).p2_u;u3=p2_data(3).p2_u;u4=p2_data(4).p2_u;u5=p2_data(5).p2_u;
y1=p2_data(1).p2_y;y2=p2_data(2).p2_y;y3=p2_data(3).p2_y;y4=p2_data(4).p2_y;y5=p2_data(5).p2_y;
u6=p2_data(6).p2_u;u7=p2_data(7).p2_u;u8=p2_data(8).p2_u;u9=p2_data(9).p2_u;u10=p2_data(10).p2_u;
y6=p2_data(6).p2_y;y7=p2_data(7).p2_y;y8=p2_data(8).p2_y;y9=p2_data(9).p2_y;y10=p2_data(10).p2_y;
u11=p2_data(11).p2_u;u12=p2_data(12).p2_u;u13=p2_data(13).p2_u;u14=p2_data(14).p2_u;u15=p2_data(15).p2_u;
y11=p2_data(11).p2_y;y12=p2_data(12).p2_y;y13=p2_data(13).p2_y;y14=p2_data(14).p2_y;y15=p2_data(15).p2_y;

Ts = 0.1;

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder, refleting your legi number.

% Use the vectors p2_u, p2_y to solve the problem. 

% Modify your code in the next section, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code.

% Do not print the values of variables (use ";" at the end of each line)
% Only show the plots required in the script.

% Check your solution with the function HS2022_SysID_midterm_p2_check
% before submitting.

%% Write your solutions here
%% 1
disp("===== 1 =====");
disp("(a) Pseudo Random Binary Signal. Because the PRBS is periodic and we can employ this perodicity to deal with trasient error.")
disp("Specifically, I choose ID=14. Because M=255 and to calculate FFT employing periodic signal,")
disp("we require the measurement length to be integer multiples of M. ID=14 satisfies this and has")
disp("highest signal-to-noise ratio to be robust against noise at high frequency")
disp(newline)
disp("(b) I first get the discrete fourier transform of the input signal U(e^jw) and output signal Y(e^jw).")
disp("Then performing ETFE, p21_G_hat = Y(e^jw)/U(e^jw). Noting that it's an element-wise division and")
disp("we only take the first k elements of the vector.")
disp(newline)
disp("(c)Yes, unbiased. Because p21_G_hat=Y(e^jw)/U(e^jw)=G_true(e^jw)+R(e^jw)/U(e^jw)+V(e^jw)/U(e^jw).")
disp("Here, R(e^jw) is the FFT of the trasient which can be eliminated because input is periodic.")
disp("This is because when calculating FFT of input, we know u(k) at k<0 due to periodicity.")
disp("Meanwhile, the noise is zero mean, and hence the expectation of its FFT is also zero. Therefore, the expectation of p21_G_hat is equal to G_true.")
p21_ID = [14]; %An integer in {1, ..., 15}
u21 = p2_data(p21_ID).p2_u;
y21 = p2_data(p21_ID).p2_y;
M=255;k=127;
% %plot the input and output signal
% figure(1)
% t_vec=0:Ts:Ts*(length(u21)-1);
% plot(t_vec,u21); hold on; plot(t_vec,y21); grid on; legend('u','y')

%get the length of each experiment
L=[];
for i= 1:15
    L=[L,length(p2_data(i).p2_u)];
end
%Average the input and output signal and discard the first period
y21_1=reshape(y21(M+1:end),M,[]);
y21_final=mean(y21_1,2);
u21_1=reshape(u21(M+1:end),M,[]);
u21_final=mean(u21_1,2);
Y21=fft(y21_final); U21=fft(u21_final);  %only positive frequencies

p21_omega = 2*pi/(Ts*M)*(1:k);
p21_G_hat = Y21(2:k+1)./U21(2:k+1); % throw away the zero frequency
figure(211);
subplot(2,1,1)
semilogx(p21_omega,20*log10(abs(p21_G_hat)))
title('Magnitude')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
xlim([min(p21_omega),max(p21_omega)])
subplot(2,1,2)
semilogx(p21_omega,rad2deg(angle(p21_G_hat)))
title('Phase')
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
xlim([min(p21_omega),max(p21_omega)])
%% 2
disp(newline)
disp("===== 2 =====");
disp("(a) I choose ID=7 because we require measurement length to be interger multiple of M")
disp("Both ID=7 and ID=9 satisfy this, and ID=7 has much higher input magnitude.")
disp(newline)
disp("(b) The cross spectrum phi_yu(e^jw)=G(e^jw)*phi_u(e^jw)+phi_uv(e^jw). If we use open-loop, the input and noise are uncorelated and hence phi_uv(e^jw)=0.")
disp("So, p22_G_hat=phi_hat_yu(e^jw)/phi_hat_u(e^jw). Here, both power spectrums are only estimated value.")
disp(newline)
disp("(c) No, it is biased. Because for random signal, the autocorrelation estimate is only asymptotically unbiased.")
disp("and since the number of measurements is finite, the autocorrelation estimate is biased, E{R(tau)} = (N-|tau|)/N*R_true(tau).")
disp("Because the power spectrum is calculted by taking FFT of R_yu(tau) and R_u(tau), the final estimate of G_hat2(e^jw) is biased.")
disp(newline)
disp("(d) For random input signal, its FFT can be very low at certain frequency and hence V(e^jw)/U(e^jw)")
disp("could be high at high frequency, which cause noisy oscillation of the estimate using ETFE")
disp("On the other hand, spectral method uses division of power spectrum, and RBS has a flat spectrum on average")
disp("So, if we average different experiments' input, the plant estimate result would be less oscillated.")
disp("The second reason is that EFTE may have large trasient error using RBS because it is not periodic.")
p22_ID = [7]; %An integer in {1, ..., 15}
u22 = p2_data(p22_ID).p2_u;
y22 = p2_data(p22_ID).p2_y;
y22_1=reshape(y22(M+1:end),M,[]);
y22_final=mean(y22_1,2);
u22_1=reshape(u22(M+1:end),M,[]);
u22_final=mean(u22_1,2);
N=length(u22_final);

%Autocorelation
tau=-(N-1)/2:(N-1)/2;
R_u=zeros(N,1); R_yu=zeros(N,1); 
for t=1:length(tau)
    for lamda=1:N
        idx=lamda-tau(t);
        if idx<1
            idx=lamda-tau(t)+N;
        elseif idx>N
            idx=lamda-tau(t)-N;
        end
        R_u(t,1)=R_u(t,1)+1/N*u22_final(lamda)*u22_final(idx);
        R_yu(t,1)=R_yu(t,1)+1/N*y22_final(lamda)*u22_final(idx);
    end
end

split = find(tau==0);

% %change to tau in [0:N-1]
% R_u_final=zeros(N,1);R_yu_final=zeros(N,1);
% R_u_final(1)=R_u(split);R_yu_final(split)=R_yu(split);
% for i= (split+1):N
%    R_u_final(i-split+1)=R_u(i); R_u_final(i)=R_u(i-split);
%    R_yu_final(i-split+1)=R_yu(i);R_yu_final(i)=R_yu(i-split);
% end
p22_phi_u = DFT(R_u, M); % DFT from -pi to pi
p22_phi_yu = DFT(R_yu, M);

p22_omega = 2*pi/(Ts*M)*(1:k);
p22_G_hat = p22_phi_yu(split+1:end)./p22_phi_u(split+1:end); %find positive frequency index
% X1 = sprintf('Length of R_u %d:',length(R_u));
% disp(X1)
% X2 = sprintf('Length of phi_u %d:',length(p22_phi_u));
% disp(X2)
% U22=fft(u22_final);
% figure(11)
% plot(tau,R_u); hold on; plot(tau,R_yu)
% figure(12)
% plot(tau,R_u_final); hold on; plot(tau,R_yu_final)
% figure(6)
% plot(p22_omega,abs(U22(1:k)))
% title('|U(e^{jw})| (linear scale)')
% xlabel('Frequency (rad/s)'); ylabel('Magnitude');
figure(221); 
plot(p22_omega,20*log10(abs(p22_phi_u(split+1:end))))
title('Magnitude of input PSD (linear scale)')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
figure(222);
plot(p22_omega,20*log10(abs(p22_phi_yu(split+1:end))))
title('Magnitude of cross PSD (linear scale)')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
figure(223);
subplot(2,1,1)
semilogx(p22_omega,20*log10(abs(p22_G_hat)))
title('Magnitude of plant estimate')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
xlim([min(p22_omega),max(p22_omega)])
subplot(2,1,2)
semilogx(p22_omega,rad2deg(angle(p22_G_hat)))
title('Phase of plant estimate')
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
xlim([min(p22_omega),max(p22_omega)])

% U22=fft(u22_final);
% figure(100)
% plot(p22_omega,abs(U22(1:k)))
% title('|U(e^{jw})| (linear scale)')
% xlabel('Frequency (rad/s)'); ylabel('Magnitude');
%% 3 
% [1.12,5.37]--> k1=36; k2=175;
disp(newline)
disp("===== 3 =====");
disp("(a) I use Chirp signal, specifically ID=1. Chirp signal is suitable when we specify the frequencies")
disp("within a certain interval. It provides very high resolution and high energy concentration within this interval.")
disp("Since I take FFT and M=2047, the measurement length should be interger multiple of M.")
disp("I plot the magnitude of FFT of input signal within [1.12,5.37]rad/s and find that only ID=1")
disp("has high input magnitude at all frequencies within this range, while the others have very low magnitude at frequency range [4,5.37]rad/s.")
disp("This would cause a zero-division using ETFE and the consequent estimate plot will have big oscillation at [4,5.37]rad/s")

M3=2047; k1=36; k2=175;
p23_ID = [1]; %An integer in {1, ..., 15}
u23 = p2_data(p23_ID).p2_u;
y23 = p2_data(p23_ID).p2_y;
y23_1=reshape(y23(M3+1:end),M3,[]);
y23_final=mean(y23_1,2);
u23_1=reshape(u23(M3+1:end),M3,[]);
u23_final=mean(u23_1,2);
Y23=fft(y23_final); U23=fft(u23_final);  %only positive frequencies
p23_omega = 2*pi/(Ts*M3)*(36:175);
p23_G_hat = Y23(k1+1:k2+1)./U23(k1+1:k2+1);
% figure(8)
% plot(p23_omega,abs(U23(k1:k2)))
% title('|U(e^{jw})| (linear scale)')
% xlabel('Frequency (rad/s)'); ylabel('Magnitude');
% This figure not required by the question
figure(3);
subplot(2,1,1)
semilogx(p23_omega,20*log10(abs(p23_G_hat)))
title('Magnitude')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
xlim([min(p23_omega),max(p23_omega)])
subplot(2,1,2)
semilogx(p23_omega,rad2deg(angle(p23_G_hat)))
title('Phase')
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
xlim([min(p23_omega),max(p23_omega)])
sgtitle("ETFE for frequency within [1.12,5.37]rad/s")
% Write the DFT function so that frequency from (-pi.pi)
function U=DFT(u,M)
    Nx=length(u);
    omega=2*pi/M*(-(M-1)/2:(M-1)/2);
    U=zeros(M,1);
    for n=1:M
        omega_n=omega(n);
        for k_x=1:Nx
            U(n)=U(n)+u(k_x)*exp(-j*omega_n*k_x);
        end
    end
end
end


