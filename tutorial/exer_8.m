%% Q1
% Samples
N=1e4; Ts=1;
%System
A=tf([1 -1.5 0.7],1.0,Ts,'Variable','z^-1');
B=tf([0 1.0 0.5],1.0,Ts,'Variable','z^-1');
C=tf([1 -1.0 0.2],1.0,Ts,'Variable','z^-1');
%Y=G*u+H*e;
G=B/A;  H=C/A;
% Input
L=tf([0 1 0.2],[1 -0.1 -0.12],Ts,'Variable','z^-1');

    e=randn(N,1);
    e_u=randn(N,1);
    %response
    u=lsim(L,e_u); % u=L*e_u
    y=lsim(G,u)+lsim(H,e); %Y=G*u+H*e;
    
    u_F=lsim(1/C,u);
    y_F=lsim(1/C,y);
    Phi_1=[0 0 0 0;-y_F(1) 0 u_F(1) 0;]; %assuem system at rest at t<=0
    for k=3:N
        Phi_1=[Phi_1;[-y_F(k-1), -y_F(k-2), u_F(k-1), u_F(k-2)]];
    end
    theta_1=Phi_1\y_F;
%% Q2 with validation dataset
N_est=N-1000;
 Phi_est=[0 0 0 0;-y_F(1) 0 u_F(1) 0;]; %assuem system at rest at t<=0
    for k=3:N_est
        Phi_est=[Phi_est;[-y_F(k-1), -y_F(k-2), u_F(k-1), u_F(k-2)]];
    end
    theta_est=Phi_est\y_F(1:N_est);
    
%% Q3 repeat 100 realizations
%Estimation
TH=[];
for l=1:100 % 100 different realizations
    e=randn(N,1);
    e_u=randn(N,1);
    %response
    u=lsim(L,e_u); % u=L*e_u
    y=lsim(G,u)+lsim(H,e); %Y=G*u+H*e;
    
    u_F=lsim(1/C,u);
    y_F=lsim(1/C,y);
    Phi=[0 0 0 0;-y_F(1) 0 u_F(1) 0;]; %assuem system at rest at t<=0
    for k=3:N_est(end)
        Phi=[Phi;[-y_F(k-1), -y_F(k-2), u_F(k-1), u_F(k-2)]];
    end
    theta=Phi\y_F(N_est);
    TH=[TH,theta];
end
%%
A=tf([1 -1.5 0.7],1.0,'Variable','z^-1');
B=tf([0 1.0 0.5],1.0,'Variable','z^-1');
C=tf([1 -1.0 0.2],1.0,'Variable','z^-1');
%Y=G*u+H*e;
G=B/A;  H=C/A;