
clc; clear all; close all

[p21_ab, p21_tau_max, p22_ghat, p22_covg, p23_ghat,p2_u,p2_y,Phi_u] = HS2022_SysID_final_p2_21957907();
%%
A=tf([1 a1 a2],1.0,Ts,'Variable','z^-1');
B=tf([0 b1],1.0,Ts,'Variable','z^-1');
G_hat=B/A;
y_hat=lsim(G_hat,p2_u);
err=p2_y-y_hat;
figure(1)
plot([1:length(p2_y)],p2_y)
hold on
plot([1:length(p2_y)],y_hat)
legend('p2_y','y_{hat}')
figure(2)
plot([1:length(p2_y)],err)
%%
Ts=1;
a1=p21_ab(1); a2=p21_ab(2); b1=p21_ab(3);
A=tf([1 a1 a2],1.0,Ts,'Variable','z^-1');
B=tf([0 b1],1.0,Ts,'Variable','z^-1');
G_hat=B/A;
t_max=0:100;
g_hat=impulse(G_hat,t_max);
N=length(p2_u); n_g=length(g_hat);
E=zeros(N,1);
% tau=13; k=38;
% e=trun_err(k,tau,n_g,g_hat,p2_u)
% for i=1:N
%  E(i)=trun_err(i,tau,n_g,g_hat,p2_u);
% end 
% e_max=max(E)

for tau=1:100
    E=zeros(N,1);
    for i=1:N
        E(i)=trun_err(i,tau,n_g,g_hat,p2_u);
    end 
    e_max=max(E);
    if e_max<=0.01
       disp(tau);
       disp(e_max)
       break; 
    end
end
% 从index 1开始计算
function e=trun_err(k,tau,n_g,g_hat,u) % n_g is the length of g_hat (pulse response)
e=0;
    if k<=tau;
        e=0;
    elseif k < n_g 
        e=0;
         for i=(tau+1):k
            e=e+g_hat(i)*u(k-i+1);
        end
    else 
        e=0;
         for i=(tau+1):n_g
            e=e+g_hat(i)*u(k-i+1);
        end
    end
end
