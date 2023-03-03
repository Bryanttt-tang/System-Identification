function [p3_Tur, p3_ID, p3_beta, p3_alpha, p3_Rhat, p3_Ghat] = HS2022_SysID_final_p3_21957907()
    %% Solution template for Problem 3
    
    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder.
    
    % Modify your code in the next sections, and return the variables
    % requested.
    
    % If you skip one part of the problem, return the null variables as already
    % provided in the code.
    
    % Please do NOT use "Run Section" to run the script.
    
    % Make sure to run the check function before submiting.
    
    % Extract Legi from Filename
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    % Obtain experiment data
    [p3_wl,p3_wu,p3_a,p3_r,p3_y] = HS2022_SysID_final_p3_GenerateData(LegiNumber);
    
    N = size(p3_r,1);
    tau_max = 20;
    
    %% Task 1
    disp("===== 1 =====");
    disp("(a)Since we are only given data of p3_r and p3_y, we first find the ETFE estimate of Tyr=Y(e^jw)./R(e^jw)")
    disp("To perform ETFE, we first determine the period of three reference signals using auto-correlation. We find that all three signals have period=500 with data length also equal to 500 ")
    disp("So, the number of period is 1 and we cannot remove the first period of r(k) and y(k) because otherwise we won't have any data left.")
    disp("Because we assume the system was at rest, therefore under this assumption, the trasients can be neglected")
    disp("Then we do element-wise division of the FFTs of y(k) and r(k) at the required frequencies to get Tyr=Y(e^jw)./R(e^jw")
    disp("Because the relation between u(k) and r(k) is: u(k)=r(k)-Co(z)*y(k). Thus,for each experiment, Tur=U(e^jw)./R(e^jw)=1-Co.*Tyr")
    disp(", where Co is the frequency response of the controller Co(z) at the required frequencies and it is a pointwise multiplication of Co and Tyr")
    disp("Alternatively, because u(k)=r(k)-Co(z)*y(k), and Co(z) is given, we can use lsim to get the time-domain data of u(k).")
    disp("And then similar to the above, we do FFTs of u(k) and r(k) to get Tur=U(e^jw)./R(e^jw)")
    disp(newline)
    disp("(b)Yes, it is asymptotically unbiased. Because Tur=1-Co.*Tyr, where Co is a known controller and can be treated as an exactly known frequency response.")
    disp("Thus, Tur is linear in Tyr. Tur is asymptotically unbiased if and only if Tyr is asymptotically unbiased.")
    disp("Note that the noise is zero mean and uncorrelated with r(k) (i.e.,Phi_er=Phi_vr=0, cross spectrum of noise and r(k) is zero), thus Tyr=Y(e^jw)./R(e^jw) is asymptotically unbiased, so as is p3_Tur.")
    disp("If we use the other method: directly get u(k) via lsim and then perform ETFE to get Tur=U(e^jw)./R(e^jw), the estimate is still asymptotically unbiased because the noise is zero mean and uncorrelated with r(k). ")
    disp("Note that the G_hat estimate of the plant is G_hat=Tyr./Tur is biased because the noise enters in a nonlinear manner via the spectra of S and S*Co.")
    disp(newline)
    disp("(d) We first compare the sensitivity function S of three experiments at the interested frequency [wl,wu]. The sensitivity function can be obtained by observing the Tur.")
    disp("I found that both experiment 1 and 2 have fairly large magnitude of S at the range [wl,wu], whereas experiment 3 has much lower magnitude.")
    disp("Large magnitude of sensitivity function means more information about the system at [wl,wu] and hence better identification performance at this frequency range. So experiment 3 is discarded.")
    disp("Then, I compare the power spectrum density of the excitation signal r of experiment 1 and 2 at [wl,wu]. Because higher PSD means higher signal-to-noise ratio, which leads to better estimation.")
    disp("I found that experiment 2 has larger PSD at the given range. So, experiment 2 would obtain the smallest mean square error.")
    
    Ts=1; k=N/2;
    a1=p3_a(1); a2=p3_a(2);a3=p3_a(3);
    C_0z_1 = tf(a1, [1 -1], Ts); C_0z_2 = tf(a2, [1 -1], Ts); C_0z_3 = tf(a3, [1 -1], Ts);

Y1=fft(p3_y(:,1)); Y2=fft(p3_y(:,2)); Y3=fft(p3_y(:,3));
R1=fft(p3_r(:,1)); R2=fft(p3_r(:,2)); R3=fft(p3_r(:,3)); %only positive frequencies

omega_k= 2*pi/(Ts*N)*(1:k);
  p3_Tur      = zeros(N/2,3);             % matrix of dimension N/2 x 3
%method 1
T_yr1 = Y1(2:k+1)./R1(2:k+1); T_yr2 = Y2(2:k+1)./R2(2:k+1); T_yr3 = Y3(2:k+1)./R3(2:k+1);
C_fre_1=squeeze(freqresp(C_0z_1,omega_k)); C_fre_2=squeeze(freqresp(C_0z_2,omega_k)); C_fre_3=squeeze(freqresp(C_0z_3,omega_k));
% p3_Tur(:,1)=1-C_fre_1.*T_yr1; p3_Tur(:,2)=1-C_fre_2.*T_yr2; p3_Tur(:,3)=1-C_fre_3.*T_yr3;
%method 2 (being used to get p3_Tur)
r1=p3_r(:,1); r2=p3_r(:,2); r3=p3_r(:,3);
u1=p3_r(:,1)-lsim(C_0z_1,p3_y(:,1)); u2=p3_r(:,2)-lsim(C_0z_2,p3_y(:,2)); u3=p3_r(:,3)-lsim(C_0z_3,p3_y(:,3));
U1=fft(u1); U2=fft(u2); U3=fft(u3);
% sensitivity function at [wl,wu]
S1=U1(2:k+1)./R1(2:k+1); S2=U2(2:k+1)./R2(2:k+1); S3=U3(2:k+1)./R3(2:k+1);
p3_Tur(:,1)=S1; p3_Tur(:,2)=S2; p3_Tur(:,3)=S3;

T_ur1=p3_Tur(:,1); T_ur2=p3_Tur(:,2); T_ur3=p3_Tur(:,3);
G_1=T_yr1./T_ur1; G_2=T_yr2./T_ur2; G_3=T_yr3./T_ur3;

omega_full= 2*pi/(Ts*N)*(1:N);

figure(31)
semilogx(omega_k,20*log10(abs(p3_Tur(:,1)))); hold on;
semilogx(omega_k,20*log10(abs(p3_Tur(:,2)))); hold on;
semilogx(omega_k,20*log10(abs(p3_Tur(:,3)))); hold off; grid on;
title('p3_{Tur}')
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
xlim([min(omega_k),max(omega_k)])
legend('Experiment 1','Experiment 2','Experiment 3')

% Place at [wl,wu]
range=find(omega_k>=p3_wl & omega_k<=p3_wu);
range_plus=[range,max(range)+1]; % include wu

% sensitivity function at [wl,wu]
% figure(100)
% semilogx(omega_k(range_plus),20*log10(abs(S1(range_plus)))); hold on;
% semilogx(omega_k(range_plus),20*log10(abs(S2(range_plus)))); hold on;
% semilogx(omega_k(range_plus),20*log10(abs(S3(range_plus)))); hold off; grid on;
% title('S')
% xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
% % xlim([min(omega_k),max(omega_k)])
% legend('Experiment 1','Experiment 2','Experiment 3')
% R11=fft(r1); R22=fft(r2); R33=fft(r3);
% p1=abs(R11(2:k+1)).^2/N; p2=abs(R22(2:k+1)).^2/N; p3=abs(R33(2:k+1)).^2/N;
% % PSD of excitation r at the given range
% figure(250) 
% plot(omega_k(range_plus),p1(range_plus)); hold on;
% plot(omega_k(range_plus),p2(range_plus)); hold on;
% plot(omega_k(range_plus),p3(range_plus)); hold off; grid on
% title('PSD')
% legend('Experiment 1','Experiment 2','Experiment 3')

% %Autocorelation
% tau_auto=-N/2+1:N/2;
% R_r_1=zeros(N,1); R_r_2=zeros(N,1); R_r_3=zeros(N,1);
% for t=1:length(tau_auto)
%     for lamda=1:N
%         idx=lamda-tau_auto(t);
%         if idx<1
%             idx=lamda-tau_auto(t)+N;
%         elseif idx>N
%             idx=lamda-tau_auto(t)-N;
%         end
%         R_r_1(t,1)=R_r_1(t,1)+1/N*r1(lamda)*r1(idx);
%         R_r_2(t,1)=R_r_2(t,1)+1/N*r2(lamda)*r2(idx);
%         R_r_3(t,1)=R_r_3(t,1)+1/N*r3(lamda)*r3(idx);
%     end
% end
% split = find(tau_auto==0);
% phi_r1= DFT(R_r_1, N); phi_r2= DFT(R_r_2, N); phi_r3= DFT(R_r_3, N);
% 
% figure(251)
% plot(omega_k(range_plus),phi_r1(split+range_plus)); hold on;
% plot(omega_k(range_plus),phi_r2(split+range_plus)); hold on;
% plot(omega_k(range_plus),phi_r3(split+range_plus)); hold off; grid on
% title('PSD')
% legend('Experiment 1','Experiment 2','Experiment 3')


    p3_ID       = 2;                        % scalar in {1,2,3}
    
    
    %% Task 2
    disp("===== 2 =====");
    disp("(a)By Youla parametrisations, we know that for a given controller C=Xo/Yo, all plants stable in closed-loop with this controller has the form:")
    disp("G_R=N/D=(No+R*Yo)/(Do-R*Xo), where R is stable. Meanhwhile, the noise model H(z) can be modeled by H=F/D=F/(Do-R*Xo), where F is stable and stably invertible.")
    disp("Then, the closed-loop system is equivalent to an open-loop coprime model y=N/D*u+F/D*e. Doing substitution, we get")
    disp("(Do-R*Xo)*y = (No+R*Yo)*u+F*e. By rearranging the terms, Do*y-No*u=R*(Xo*y+Yo*u)+F*e. So, beta=Do*y-No*u and alpha=Xo*y+Yo*u, with R and F stable.")
    disp("Since we are only given p3_y and p3_r, we have to do further steps. With u=r-Co*y, we can do substitution")
    disp("Eventually, we get beta=(Do+No*Xo/Yo)*y-No*r and alpha=Yo*r. Then we can use lsim command to construct p3_beta and p3_alpha")
    disp(newline)
    disp("(b) No, it is biased. Because to estimate G_hat, we need to estimate R_hat and let G_hat=(No+R_hat*Yo)/(Do-R_hat*Xo)")
    disp("Here, R_hat estimate is already biased because we estimate its pulse response using finite pulse response estimate.")
    disp("It will have truncation error and will lead to biased estimate of R_hat, which inevitably causes G_hat biased.")
    disp("Furthermore, even if the R_hat is unbiased, the ratio (No+R_hat*Yo)/(Do-R_hat*Xo) can also be biased because it is a nonlinear operation of two unbiased estimates, which generally results in biased estimate")
    disp(newline)
    disp("(c) i. False. G_hat is stablized by the controller. Because by Dual Youla method, any G_hat of the form:")
    disp("G_hat=(No+R_hat*Yo)/(Do-R_hat*Xo) is stablized by C0 as long as R_hat is stable and the initial plant P0=N0/D0 is stable in closed-loop with C0.")
    disp("I construct the initial closed loop transfer function and test out that it is indeed stable with controller C0.")
    disp("Meanwhile, the finite pulse response estimate R_hat is stable. In fact, R is stable because only if R is stable can we find tau_max to estimate FIR. Thus, G_hat is stable")
    disp(newline)
    disp("ii. False. R_hat is stable because it is FIR estimate. We can verify this emperically by plotting the response of R_hat and finding that it converges to zero.")
    disp("In fact, R is stable because only if R is stable can we find tau_max to estimate FIR. Thus, a FIR estimate of a stable R is also stable.")
    disp(newline)
    disp("iii. True. This is the same reasoning as question (i). Because G_hat=(No+R_hat*Yo)/(Do-R_hat*Xo), G_hat is stable if R_hat is stable.")
    
    N0=tf(1-0.2, [1 -0.2], Ts); D0=1; P0=N0/D0;
     a=p3_a(p3_ID); r=p3_r(:,p3_ID); y=p3_y(:,p3_ID);
     C_0 = tf(a, [1 -1], Ts);
     sys_initial=C_0*P0/(1+C_0*P0);
%      figure(111)
%      step(sys_initial); % test if P0 is initally stabilized by C0
     [fact,Y0,X0] = rncf(C_0);
     beta1=D0+N0*X0/Y0;
    p3_beta     = lsim(beta1,y)-lsim(N0,r); % vector of dimension N x 1
    p3_alpha    = lsim(Y0,r);               % vector of dimension N x 1
    
    Phi_u=zeros(N,tau_max+1);
for i=1:N
   if i<=tau_max+1
       for j=1:i
        Phi_u(i,j)=p3_alpha(i-j+1); 
       end
   else
       for j=1:(tau_max+1)
        Phi_u(i,j)=p3_alpha(i-j+1);
       end
   end    
end
    p3_Rhat     = Phi_u\p3_beta;       % vector of dimension tau_max+1 x 1
R=0; 
z=tf('z',1);
for i=1:length(p3_Rhat)
    R=R+p3_Rhat(i)*z^(-i+1);
end
  G=(N0+R*Y0)/(D0-R*X0);
  
    p3_Ghat = squeeze(freqresp(G,omega_k));             % vector of dimension N/2 x 1
% figure(1000)
% semilogx(omega_k,20*log10(abs(G_1))); hold on;
% semilogx(omega_k,20*log10(abs(G_2))); hold on;
% semilogx(omega_k,20*log10(abs(G_3))); hold on;
% semilogx(omega_k,20*log10(abs(p3_Ghat))); hold on;grid on;
% title('G')
% xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
% % xlim([min(omega_k),max(omega_k)])
% legend('Experiment 1','Experiment 2','Experiment 3','Youla')
    
end

