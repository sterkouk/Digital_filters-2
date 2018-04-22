% C 3o Erwthma
clear all
close all
tic     
% The same steps of initializations and creations
% as the previous two exercises
M = 1024; 
L = M;     
s2d = 0.57; 
n = 1.2e6; 
kmax = floor(n/L); 
u = zeros(n,1);
W = zeros(2*M,1); % Extra matrix
P = zeros(2*M,1); % Extra matrix
u_help = zeros(2*M,1); % Extra matrix
y_help = zeros(n,1); % Extra matrix
% Those were created in order to help me calculate the real ones
e = zeros(n,1);
J = zeros(n,1);
average_J = zeros(kmax-1);
v = sqrt(s2d)*randn(n,1); 
v = v - mean(v);
u(1) = v(1);
for i=2:n
    u(i) = (-0.34*u(i-1)) + v(i);
end  
d = plant(u')';

for k=1:kmax-1
    u_help = u((k-1)*M+1:(k+1)*M);
    up = fft(u_help, 2*M);
    y_help = ifft(up.*W);
    y = y_help(M+1:2*M);
    dp = d(k*M+1:(k+1)*M);
    e(k*M+1:(k+1)*M ,1) = dp-y;
    average_J(k) = sum(e(k*M+1:(k+1)*M,1).^2)/M;
    % Calculations for FFT based on the lessons the notes 
    % and an internet search
    % http://www.eit.lth.se/fileadmin/eit/courses/ett042/CE/CE4e.pdf
    U = fft([zeros(M,1)' e(k*M+1:(k+1)*M)']' ,2*M );
    g = 0.3;
    P = g*P+(1-g)*abs(up).^2;
    D = 1./P;           
    a_help = ifft(D.*conj(up).*U,2*M);
    a = a_help(1:M);
    step = 0.4; 
    W = W+step*fft([a;zeros(M,1)],2*M);
end
wf = ifft(W);
wf = real(wf(1:length(W)/2));

toc

figure(3)
semilogy(average_J);
xlabel('k')
ylabel('Ee^{2}(n)');
title('Learning curves 3rd Diagram');
