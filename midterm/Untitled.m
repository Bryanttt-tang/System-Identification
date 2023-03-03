N=6;  
sigma_e=zeros(N,N);
    for i=1:N
        sigma_e(i,i)=0.24; sigma_e(i,i+1)=0.05;sigma_e(i+1,i)=0.05;
        sigma_e(i,i+2)=0.15;sigma_e(i+2,i)=0.15;
    end