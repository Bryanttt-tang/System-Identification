function [p21_ab, p21_tau_max, p22_ghat, p22_covg, p23_ghat] = HS2022_SysID_final_p2_21957907()

%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[p2_u, p2_y] = HS2022_SysID_final_p2_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder, refleting your legi number.

% Use the vectors p2_u, p2_y to solve the problem. 

% Modify your code in the next section, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code.

% Do not print the values of variables (use ";" at the end of each line)
% Only show the plots required in the script, if any.

% Check your solution with the check function before submitting.

%% Write your solutions here
%% 1
disp("===== 1 =====");
disp("(a) This is the ARX (equation error model), A(Z)*y(k)=B(z)*u(k)+e(k). Here, A(z)=1+a1*z^-1+a2*z^-2, and B(z)=b1*z^-1")
disp("and the noise e(k) is assumed to be independent and identically distributed. Let G(z)=B(z)/A(z) and H(z)=1/A(z), where A(z) is monic and G(z) is strictly proper.")
disp("Therefore, the equation error can be written as: y_hat=1/H(z)*G(z)*u(k)+(1-1/H(z))*y(k)=B(z)*u(k)+(1-A(z))*y(k)=Phi'*p21_ab. This relation is linear in p21_ab.")
disp("Here, a row of the regressor Phi(k) would be [-p2_y(k-1), -p2_y(k-2), p2_u(k-1)], and regressor Phi would be [Phi(1); Phi(2); ...; Phi(N)]")
disp("Because we assume the system is at rest at k<=0, in the regressor formulation, y(-1)=y(0)=0; u(0)=0.")
disp("Because the AR dynamics are the same for both input and noise, the output y vector and the regressor have exactly the same noise and hence it is a linear regression problem.")
disp("Because the noise is assumed to be i.i.d, we can find p21_ab by a least square formulation and p21_ab=inv(Phi'*Phi)*Phi'*p2_y")
disp(newline)
disp("(b) Because this is a stable and causal system, for any epsilon >0, we can find tau_max such that summation of g(i) for i from tau_max+1 to infinity will < epsilon")
disp("From p21_ab, we can get the estimated approximate transfer function model using tf command. Then using impulse() command, we get the pulse response estimate g_hat")
disp("Then, I write a function trun_err to calculate the truncation error, which is basically a convolution operation. Then I loop for different tau_max values, I find max(abs(e(k)))<=0.01")
disp("Note that the MATLAB starts from index 1 while in the question, g starts from index 0. So, p21_tau_max=return value-1")

N=length(p2_y);
Phi=zeros(N,3);
Phi(1,:)=[0 0 0];
Phi(2,:)=[-p2_y(1) 0 p2_u(1)];
for k=3:N
    Phi(k,:)=[-p2_y(k-1), -p2_y(k-2), p2_u(k-1)];
end

p21_ab = Phi\p2_y;

Ts=1;
a1=p21_ab(1); a2=p21_ab(2); b1=p21_ab(3);
A=tf([1 a1 a2],1.0,Ts,'Variable','z^-1');
B=tf([0 b1],1.0,Ts,'Variable','z^-1');
G_hat=B/A;
t_max=0:100;
g_hat=impulse(G_hat,t_max); n_g=length(g_hat);

for tau=1:100
    E=zeros(N,1);
    for k=1:N
        E(k)=trun_err(k,tau,n_g,g_hat,p2_u);
    end 
    e_max=max(abs(E));
    if e_max<=0.01
       t_cur=tau;
       disp(e_max)
       break; 
    end
end
% p21_tau_max = 14;
p21_tau_max=t_cur-1;

%% 2
disp("===== 2 =====");
disp("(a) For the finite pulse response estimate, we can only estimate tau_max+1 parameters of g_hat.")
disp("From the pulse response model, we can write the convolution of g and u into matrix formulation: p2_y=Phi_u* p22_ghat+V, where V is a column of zero-mean i.i.d noise.")
disp("Here, the regressor Phi_u is a N-by-tau_max+1 matrix, where N is the length of output measurement p2_y. Each row of the regressor Phi_u(k)=[u(k),u(k-1),..,u(k-tau_max)].")
disp("So, Phi_u=[Phi_u(1); Phi_u(2); ...; Phi_u(N)]. Noting that since N> tau_max+1, the regressor is a tall matrix.")
disp("Because we assume the system is at rest at k<=0, in the regressor formulation, u(k)=0 for k<=0.")
disp("So, the formulation is p2_y=Phi_u* p22_ghat+V, which is linear in p22_ghat and can be treated as least square problem.")
disp("The least square solution is p22_ghat=inv(Phi_u'*Phi_u)*Phi_u'*p2_y. Note that we need (Phi_u'*Phi_u) to be invertible (i.e.,u is persistently exciting)")
disp("Note that because the noise is zero-mean, the p22_ghat is unbiased as long as there is no trasient error.")
disp(newline)
disp("(b) If we consider the truncation error, the whole formulation is p2_y=Phi_u* p22_ghat+V+e where e is the trunction error vector.")
disp("Here, e is not a stochastic vector and will only contribute to the bias; whereas V is a zero-mean stochastic vector and will contribute to the variance")
disp("Because v(k) is zero-mean and i.i.d., cov(g)=E{(g_hat-E{g_hat})*(g_hat-E{g_hat})'}=sigma^2*inv(Phi_u'*Phi_u), where sigma is the standard deviation of v(k)")

%Regressor
Phi_u=zeros(N,p21_tau_max+1);

for i=1:N
   if i<=p21_tau_max+1
       for j=1:i
        Phi_u(i,j)=p2_u(i-j+1); 
       end
   else
       for j=1:(p21_tau_max+1)
        Phi_u(i,j)=p2_u(i-j+1);
       end
   end
    
end

p22_ghat = Phi_u\p2_y;
sigma=0.25; % noise v(k)
p22_covg = sigma^2*inv(Phi_u'*Phi_u);

%% 3
disp("===== 3 =====");
disp("(a) Because the measurements with value equal to zero are corrupted measurements, we need to find the place of these values in p2_y using command find(p2_y==0)")
disp("Then we need  to delete these corrputed measurements in p2_y. Meanwhile, in the regressor matrix Phi_u, the corresponding rows k at which p2_y(k)=0 need to be deleted.")
disp("Thereby, let N be the length(p2_y) and n be the number of corrupted measurements, then we effectively use (N-n) measurements.")
disp("The new regressor matrix becomes a (N-n)-by-tau_max+1 matrix and we can still estimate p23_ghat as long as the remaining (N-n)-by-tau_max+1 matrix is full-column rank")
disp("Noting that since only the measurements are corrupted while the inputs are fine, we just need to eliminate the rows in the regressor corresponding to the place of bad measurements")

bad=find(p2_y==0); % place of corrupted measurements

Phi_u_3=Phi_u;
Phi_u_3(bad,:)=[];
p2_y_3=p2_y;
p2_y_3(bad,:)=[];

p23_ghat = Phi_u_3\p2_y_3;

%%
% index from 1. Means the first tau terms of g.
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
end


