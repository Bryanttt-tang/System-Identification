mu=0.2; sigma=0.1; mu_th=[1;0.4]; sigma_th=0.01;
W=[x,2*x.^2-1];
theta_ml= (W'*W)\W'*(y-mu);
theta_map=(1/sigma*W'*W+1/sigma_th*eye(2))\(1/sigma*W'*(y-mu)+mu_th/sigma_th);
W_val=[x_v,2*x_v.^2-1];
y_1=W_val*theta_ml+mu;
y_2=W_val*theta_map+mu;
e1=1/length(y_v)*norm(y_v-y_1)
e2=1/length(y_v)*norm(y_v-y_2)