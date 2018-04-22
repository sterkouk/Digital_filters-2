% C 1o Erwthma
clear all
close all

tic
% initialize the filter parametrs
M = 1024;
L = M;     
s2d = 0.57;
n = 1.2e6; 
mu = 0.0005;
kmax = floor(n/L);
% initialize
u = zeros(n,1);
w = zeros(M,1);
y = zeros(n,1);
e = zeros(n,1);
J = zeros(n,1);
average_J = zeros(kmax-1);
% Noise creation
v = sqrt(s2d)*randn(n,1); 
v = v-mean(v);
u(1) = v(1);
for i=2:n
    u(i) = (-0.34*u(i-1))+v(i);
end
d = plant(u')';
% Implementation of LMS algortithm
% Based on the example adaptLMS and the example from the notes
for k=1:kmax-1
    a = zeros(M,1);
    up = zeros(M,1);
    for i=1:L-1
        up = u(k*L+i:-1:k*L+i-M+1);
        dp = d(k*L+i);     
        y = w'*up;       
        e = dp-y;        
        a = a+(mu*e*up);
        J = e^2;
        average_J(k) = average_J(k)+J;
    end
    average_J(k) = average_J(k)/L;
    w = w + a;
end     

toc    

figure(1)
semilogy(average_J);
xlabel('k')
ylabel('Ee^{2}(n)');
legend({sprintf('mu=%f',mu)});
title('Learning curves 1st Diagram ');


