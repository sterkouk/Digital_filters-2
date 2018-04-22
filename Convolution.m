% B Erwthma
clear all
close all

n=6;
% Signal initialization
signal_x = [n,1];
signal_x = [zeros];
signal_y = [n,1];
signal_y = [zeros];
signal_x = randi([10 50],1,n);
signal_y =randi([10 50],1,n);
% 1o Erwthma
convolution = conv(signal_y,signal_x)
figure(1)
plot(convolution)
% 2o Erwthma
x=signal_x;
y=signal_y;
r = [y zeros(1,length(x)-1)];
c = [y(1) zeros(1,length(x)-1)];
yC = toeplitz (c,r);
apotelesma_2 = x*yC
% 3o Erwthma
x_Cirulant=signal_x;
y_Circulant=signal_y;
c_Circulant = [y_Circulant(1) fliplr(y_Circulant(2:end))];
r_Circulant = y_Circulant;
y_Circulant_conv = toeplitz(c,r);
apotelesma_3 = x_Cirulant*y_Circulant_conv
% 4o Erwthma
n1 = length(signal_x) + length(signal_y)-1 ;
fourier_x = fft(signal_x,n1);
fourier_y = fft(signal_y,n1);
apotesma_4 = ifft(fourier_x.*fourier_y)