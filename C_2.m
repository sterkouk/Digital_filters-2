% C 2o Erwthma
% Same initializations and Creations with the 1st one 
clear all
close all
tic    
M = 1024; 
L = M;      
s2d = 0.57;
n = 1.2e6;  
mu = 0.0005;
kmax = floor(n/L);
u = zeros(n,1);
w = zeros(M,1);
y = zeros(n,1);
e = zeros(n,1);
J = zeros(n,1);
average_J = zeros(kmax-1);
v = sqrt(s2d)*randn(n,1); 
v = v - mean(v);
u(1) = v(1);
for i=2:n
    u(i) = (-0.34*u(i-1))+v(i);
end
d = plant(u')';
% From here it begins the different implementation of the LMS algorithm
% Always based on the adaptLMS example and the notes
for k=1:kmax-1
    a = zeros(M,1);
    up = zeros(M,M);
    for m=0:M-1
        up(m+1,:) = u((k*M)+m:-1:((k*M)+m-M+1));
    end
    dp = d(k*M:1:(k+1)*M-1);
    y = up*w;
    e = dp-y;
    a = mu*up'*e;
    w = w+a;
    average_J(k) = sum(e.*e)/M;

end     

toc 

figure(2)
semilogy(average_J);
xlabel('k')
ylabel('Ee^{2}(n)');
legend({sprintf('mu=%f',mu)});
title('Learning curves 2nd Diagram');
