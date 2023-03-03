function [p1_theta1, p1_theta2, p1_opt, p1_y_hat] = HS2022_SysID_final_p1_21957907()
% Return vairables
% p1_theta1 - vector of dimension [4 x 1]
% p1_theta2 - vector of dimension [4 x 1]
% p1_opt    - integer between 1 and 3
% p1_y_hat  - vector of dimension [N x 1]


%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[p1_w, p1_v, p1_u, p1_y] = HS2022_SysID_final_p1_GenerateData(LegiNumber);


%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder, refleting your legi number.

% Use the vectors p1_w, p1_v, p1_u and p1_y to solve the problem. 

% Modify your code in the next section, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code.

% Do not print the values of variables (use ";" at the end of each line)
% Only show the plots required in the script.

% Check your solution with the function HS2022_SysID_final_p1_check
% before submitting.

%% Write your solutions here
%% 1
disp("===== 1 =====");
disp("(a) From the covarince between e at different time steps, we know that e(k) is correlated noise. This is the ARX (equation error model).")
disp(" We use one-step ahead prediction error based identification to obtain an estimate.")
disp("Here, the G(z)=B(z)/A(z) and H(z)=1/A(z), where A(z) is monic and G(z) is strictly proper. Therefore, the equation error can be written as")
disp("y_hat=1/H(z)*G(z)*u(k)+(1-1/H(z))*y(k)=B(z)*u(k)+(1-A(z))*y(k)=phi'*p1_theta1. This relation is linear in p1_theta1.")
disp("Here, a row of the regressor phi(k) would be [-p1_v(k-1), -p1_v(k-2), p1_w(k-1), p1_w(k-2)], and regressor Phi would be [phi(1); phi(2); ...; phi(N)]")
disp("Because we assume the system is at rest at k<=0, in the regressor formulation, w(-1)=w(0)=0; v(-1)=v(0)=0.")
disp("Because the AR dynamics are the same for both input and noise, the p1_v vector and the regressor have exactly the same noise and hence it is a linear regression problem.")
disp("Then, this problem becomes p1_v=Phi*p1_theta1+epsilon, where epsilon is the prediction error vector and epsilon(k)=e(k) for one-step ahead prediction.")
disp("Because we know noise covariance information, the BLUE estimator can be used which has the same form as the weighted least squares linear estimate (or maximum likelihood estimate).")
disp("Thus, the regression problem becomes a weighted least square problem to minimize the variance of p1_theta1 and the estimate p1_theta1=inv(Phi'*inv(sigma_e)*Phi)*Phi'*inv(sigma_e)*p1_v")
disp("Here, sigma_e is the covarince matrix of noise vector E = [e(1); e(2); ...; e(N)].")
disp("The estimate p1_theta1=Z'*p1_v, where Z'*Phi=I. Since Z' matches the expression of BLUE: Z'=(Phi'*inv(sigma_e)*Phi)\Phi'*inv(sigma_e), the estimator is asymptotically unbiased and minimum variance.")

disp("(b) p1_theta1=inv(Phi'*inv(sigma_e)*Phi)*Phi'*inv(sigma_e)*p1_v")

N1=length(p1_v);
%p1_w input
%p1_v output

% Get the covariance of E
sigma_e = zeros(N1, N1);
for i = 1:N1
    for j = 1:N1
        sigma_e(i,j) = Cov_E(i,j);
    end
end

Phi_1=zeros(N1,4);
Phi_1(1,:)=[0 0 0 0];
Phi_1(2,:)=[-p1_v(1) 0 p1_w(1) 0];
% Phi_1=[0 0 0 0;-p1_v(1) 0 p1_w(1) 0]; %assuem system at rest at t<=0
for k=3:N1
    Phi_1(k,:)=[-p1_v(k-1), -p1_v(k-2), p1_w(k-1), p1_w(k-2)];
end
 p1_theta1=(Phi_1'*inv(sigma_e)*Phi_1)\Phi_1'*inv(sigma_e)*p1_v; % vector of dimension [4 x 1]


%% 2
disp("===== 2 =====");
disp("(a) No, the instrument at time k is correlated with the noise v(k). This is because y(k)=D(z)/C(z)*u(k)+1/C(z)*v(k).")
disp("Here, y(k) is corelated with v(k), while y(k-3)is corelated with v(k-3). And from eq.(1), v(k) is corelated wtih e(k) while v(k-3) is corelated with e(k-3).")
disp("From the covariance information of e(k), we know e(k) is corelated with e(k-3)")
disp("Therefore, v(k) is corelated with v(k-3). Because y(k-3)is corelated with v(k-3), we get y(k-3) is corelated with v(k). Hence, the instrument at time k is correlated with the noise v(k).")

disp("(b)Option 3. The system has an output error structure. To obtain a consistent estimate, we need zeta(k), the instrumental variables, to be independent of the zero-mean noise epsilon(k). All three options satisfy this requirement.")
disp("Meanwhile, we need E{Zeta'*Phi} to be nonsingular, which means the instrument variables have to be sufficiently corelated with the regressor Phi.")
disp("We assume that the input u(k) is independent and identically distributed, i.e., E{u(k)*u(k+tau)}=0 for any tau not equal to 0.")
disp("From the structure of C(z) and D(z), we know that current output y(k) only depends on past input all the way up to u(k-1). The system is causal.")
disp("In other words, for k1<=k2, E{y(k1)*u(k2)}=0. Therefore, for option 1, where Zeta is N-by-4 matrix formed by stacking [Ksi(k);Ksi(k+1);...;Ksi(k+N)]")
disp("If we calculate E{Zeta' * Phi}, which returns a 4-by-4 matrix, the third column is all-zero column.") 
disp("This is because E{y(k1)*u(k2)}=0 for k1<=k2, and  E{u(k)*u(k+tau)}=0 for any tau not equal to 0. The same reasoning applies for Option 2.")
disp("Whereas for Option 3, there is no all-zero row or column, and we are able to find input such that E{Zeta' * Phi} matrix is full rank. Thus, Option 3 might obtain consistent estimate.")


N = length(p1_u);
%Regressor phi2
Phi2 = zeros(N, 4);
Phi2(1,:) = [0 0 0 0];
Phi2(2,:) = [-p1_y(1) 0 p1_u(1) 0];
for k = 3:N
    Phi2(k,:)=[-p1_y(k-1), -p1_y(k-2), p1_u(k-1), p1_u(k-2)];
end
% Zeta_2 using option 3
Zeta2 = zeros(N, 4);
Zeta2(1,:) = [0 0 0 0];
Zeta2(2,:) = [0 0 p1_u(1) 0];
Zeta2(3,:) = [0 0 p1_u(2) p1_u(1)];
Zeta2(4,:) = [p1_y(1) 0 p1_u(3) p1_u(2)];
for k = 5:N
    Zeta2(k,:)=[p1_y(k-3), p1_y(k-4), p1_u(k-1), p1_u(k-2)];
end

p1_theta2 = (Zeta2'*Phi2)\Zeta2'*p1_y; %vector of dimension [4 x 1]


%% 3
disp("===== 3 =====");



p1_opt = 3;

%Pesudolinear regressor
Ts=1;
c1=p1_theta2(1); c2=p1_theta2(2); d1=p1_theta2(3); d2=p1_theta2(4);
C=tf([1 c1 c2],1.0,Ts,'Variable','z^-1');
D=tf([0 d1 d2],1.0,Ts,'Variable','z^-1');
sys=D/C;

p1_y_hat = lsim(sys,p1_u) ; %vector of dimension [N x 1], where N = length(p1_u);


figure(131)
plot([1:N], p1_y)
hold on
plot([1:N],p1_y_hat)
hold off
xlabel('Time Step')
ylabel('Output')
legend('Measured p1_y','Estimated p1_{y_{hat}}')

%% Function
function covij = Cov_E(i,j)
    if i == j
        covij = 0.3;
    elseif abs(i-j) <= 3
        covij = 0.05;
    else
        covij = 0;
    end
end

end


