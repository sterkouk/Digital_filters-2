% C 4o Erwthma
clear all
close all
tic

% The same steps of initializations and creations
% as the previous three exercises
M = 1024;
L = M;
n = 1.4e6;
kmax = floor(n/L);
s2d = 0.57;
u = zeros(n,1);
W = zeros(2*L,1);
y = zeros(n,1);
e = zeros(n,1);
average_J = zeros(kmax);
v = sqrt(s2d)*randn(n,1); 
v = v - mean(v);
u(1) = v(1);
for i=2:n
    u(i) = (-0.34*u(i-1))+v(i);
end
d = plant(u.');
u = u.';
W = W.';
for k=2:kmax
    if k == 2
        p = var(u((k-1)*L+1:k*L))*autocorr(u((k-1)*L+1:k*L),L-1);
        mu = 0.4*(2/(L*eigs(toeplitz(p),1,'la')));
    end
    Fu = fft(u((k-2)*L+1:k*L));
    yh = ifft(Fu.*W);
    y = yh(L+1:2*L);
    e = d((k-1)*L+1:k*L)-y;
    E = fft([zeros(1,L) e]);
    W = W+mu*conj(Fu).*E;
    average_J(k) = sum(e.^2);
end       
toc 
figure(4)
semilogy(average_J)
xlabel('k')
ylabel('Ee^{2}(n)');
title('Learning curves 4th Diagram');