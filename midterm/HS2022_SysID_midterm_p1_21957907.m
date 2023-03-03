function [p1_C_ML,p1_theta_bd,p1_Sigma_bd,p1_C_MAP,p1_theta_bdMAP] = HS2022_SysID_midterm_p1_21957907()

    %% Solution template for Problem 1

    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder.

    % Modify your code in the next sections, and return the variables
    % requested.

    % If you skip one part of the problem, return the null variables as already
    % provided in the code.
    
    % Please do NOT use "Run Section" to run the script.

    % Extract Legi from Filename
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    % Obtain experiment data
    [p1_x,p1_y,p1_yf] = HS2022_SysID_midterm_p1_GenerateData(LegiNumber);
    % Length of theta parameter
    nx = size(p1_x,2);

    %% Task 1: Maximum likelihood estimate
    N=size(p1_x,1);
   % alpha=0.24; beta=0.05; gamma=0.15
    sigma_e=zeros(N,N);
    for i=1:N-2
        sigma_e(i,i)=0.24^2; sigma_e(i,i+1)=0.05^2;sigma_e(i+1,i)=0.05^2;
        sigma_e(i,i+2)=0.15^2;sigma_e(i+2,i)=0.15^2;
    end
    sigma_e(N-1,N-1)=0.24^2;sigma_e(N,N)=0.24^2;sigma_e(N,N-1)=0.05^2;sigma_e(N-1,N)=0.05^2;
     disp("===== (a) =====");
     disp("i): Given p1_x and p1_y, and let E=[e(1),...,e(N)]', we have p1_y=p1_x*C'+E. So C_ML=argmax f(p1_y|C)=argmin -log(f(p1_y|C))")
     disp("By doing derivative with respect to C and set it to zero, we get C_ML'=inv(p1_x'*inv(sigma_e)*p1_x)*p1_x'*inv(sigma_e)*p1_y")
     disp(newline)
     disp("ii): YES, it is the best linear unbiased estimator. Because we have full knowledge about the error covariance.")
     disp("The MLE is an unbiased linear estimator which satisfies C_ML=Z'*p1_y, where Z'*p1_x=I.")
     disp("Here, Z' matches the expression of BLUE: Z'=(p1_x'*inv(sigma_e)*p1_x)\p1_x'*inv(sigma_e).")
     disp("We can also prove that cov{C_ML}=inv(p1_x'*inv(sigma_e)*p1_x) <= cov{C} for any unbiased estimator C.")
     disp(newline)

    p1_C_ML = ((p1_x'*inv(sigma_e)*p1_x)\p1_x'*inv(sigma_e)*p1_y)';          % vector of dimension 1xnx


    %% Task 2: Bias & drift estimate
    disp("===== (b) =====");
    disp("i): theta_bd__ML=argmax f(p1_yf|theta_bd) = argmin -log(f(p1_yf|theta_bd).")
    disp("We reexpress the equation as p1_yf=W* theta_bd+E; where W is calculated in the code.")
    disp("By doing derivative with respect to theta_bd and set it to zero, we get:")
    disp("theta_bd_ML=inv(W'*inv(sigma_e)*W)*W'*inv(sigma_e)*p1_yf")
    disp(newline) 
    disp("ii): The MLE is an unbiased linear estimator because theta_bd_ML can be expressed as theta_bd_ML=Z'*p1_yf, where Z'*W=I.")
    disp("In this specific case, Z'=inv(W'*inv(sigma_e)*W)*W'*inv(sigma_e).")
    disp("Hence, p1_Sigma_bd=E{(theta_bd_ML-theta_true)*(theta_bd_ML-theta_true)'}=Z'*sigma_e_b*Z")
    disp(newline)
    Nf=length(p1_yf); dt=1/200; %sampling rate
     sigma_e_b=zeros(Nf,Nf); 
    for i=1:Nf-2
        sigma_e_b(i,i)=0.24^2; sigma_e_b(i,i+1)=0.05^2;sigma_e_b(i+1,i)=0.05^2;
        sigma_e_b(i,i+2)=0.15^2;sigma_e_b(i+2,i)=0.15^2;
    end
    sigma_e_b(Nf-1,Nf-1)=0.24^2;sigma_e_b(Nf,Nf)=0.24^2;sigma_e_b(Nf,Nf-1)=0.05^2;sigma_e_b(Nf-1,Nf)=0.05^2;
    % create the W matrix s.t. p1_yf=W* p1_theta_bd+E
    W=zeros(Nf,2);
    W(:,1)=ones(Nf,1);
    for i =1:Nf
        W(i,2)=dt*(i-1);
    end
    p1_theta_bd=(W'*inv(sigma_e_b)*W)\W'*inv(sigma_e_b)*(p1_yf);% vector of dimension 2x1
    p1_Sigma_bd=(W'*inv(sigma_e_b)*W)\W'*inv(sigma_e_b)*sigma_e_b*((W'*inv(sigma_e_b)*W)\W'*inv(sigma_e_b))';
 	% matrix of dimension 2x2


    %% Task 3: Maximum a posteriori estimate
disp("===== (c) =====");
disp("i):We don't have the prior of P(C) and data p1_yf is used to calcute MLE of theta_bd, which is used as the prior")
disp("Thus, if we purely use p1_yf to calculate theta_bd_MAP, the MAP will be the same of MLE results")
disp("Therefore, I stack X=[p1_x,W_new], and theta_aug=[C;theta_bd]. Here, W_new is similar to the W in part b),")
disp("except that its first dimension=N rather that Nf. Then, we get p1_y=X*theta_aug+E.")
disp("The augmented vector theta_aug_MAP= argmax f(theta_aug|p1_y).By Bayes’ rule,")
disp("f(theta_aug|p1_y) is proportional to f(p1_y|theta_aug)*f(theta_aug)")
disp("Still, we need the prior of C for the augmented vector theta_aug, and I calculated")
disp("the covariance of C_ML and assume the prior of C to be N(C_ML,C_sigma). Hence, the prior of thata_aug")
disp("can be constructed and then theta_bd_MAP= argmin -log(f(p1_y|theta_aug)-log(f(theta_aug))")
disp("= argmin 1/2*(p1_y-X*theta_aug)'*inv(Sigma_e)*(p1_y-X*theta_aug) + 1/2*(theta_aug-theta_aug_ML)'*inv(Sigma_aug)*(theta_aug-theta_aug_ML)")
disp("By doing derivative with respect to theta_bd and set it to zero, we get:")
disp("theta_aug_MAP=(X'*inv(sigma_e)*X+inv(Sigma_aug))\(inv(p1_Sigma_aug)*theta_aug_ML+X'*inv(sigma_e)*p1_y)")
disp(newline)
disp("ii): No, it won't be the same. By stack, Y_sta=[p1_y;p1_yf], and X_sta=[p1_x,W_new;zeros(Nf,nx),W]")
disp("The noise convariance sigma_e_sta=blkdiag(sigma_e,sigma_e_b) since two experiments are independent.")
disp("Thereby, Y_sta=X_sta*[C,b,d]'+E_sta. So,[C,b,d]_CL=(X_sta'*inv(sigma_e_sta)*X_sta)\X_sta'*inv(sigma_e_sta)*Y_sta.")
disp("The returned result shows that it is different from MAP. Intuitively, the experiment p1_y")
disp("contains information of both C and theta_bd, while p1_yf contains only theta_bd. So, only using")
disp("p1_y to do MAP with data from p1_yf as the prior, will have different results from combining both")
disp("p1_y and p1_yf to do MLE.")

    W_new=zeros(N,2);
    W_new(:,1)=ones(N,1);
    for i =1:N
        W_new(i,2)=dt*(i-1);
    end

X=[p1_x,W_new]; theta_aug_ml=[p1_C_ML';p1_theta_bd];
C_sigma=(p1_x'*inv(sigma_e)*p1_x)\p1_x'*inv(sigma_e)*sigma_e*((p1_x'*inv(sigma_e)*p1_x)\p1_x'*inv(sigma_e))';
Sigma_aug=blkdiag(C_sigma,p1_Sigma_bd);
theta_aug=(X'*inv(sigma_e)*X+inv(Sigma_aug))\(inv(Sigma_aug)*theta_aug_ml+X'*inv(sigma_e)*p1_y);
p1_C_MAP=theta_aug(1:nx)'; %extract the augmented vector to get C_MAP and theta_MAP
p1_theta_bdMAP=theta_aug(nx+1:end);
%     p1_C_MAP        = p1_C_ML;          % vector of wdimension 1xnx
%     p1_theta_bdMAP 	= (W'*inv(sigma_e_b)*W+inv(p1_Sigma_bd))\(inv(p1_Sigma_bd)*p1_theta_bd+W'*inv(sigma_e_b)*p1_yf);           % vector of dimension 2x1
 8.,MNBN    % ii) MLE of stack variable (not required by question)
Y_sta=[p1_y;p1_yf]; X_sta=[p1_x,W_new;zeros(Nf,nx),W];
sigma_e_sta=blkdiag(sigma_e,sigma_e_b);
ML_sta=(X_sta'*inv(sigma_e_sta)*X_sta)\X_sta'*inv(sigma_e_sta)*Y_sta;
C_ml_sta=ML_sta(1:nx)';
the_sta=ML_sta(nx+1:end);

end
