clear all;
pkg load signal

% Sampling rate
Fs = 1000;

T = 1/Fs;
M = 4000;
t = (0:M-1)*T;

%s = 20*sin((6.28*50*t));
%s = zeros(1,M);
s = 4*randn(1,M);
x = (5*sin((6.28*5*t))+5).*sin((6.28*25*t));
%x_d = (2*sin((6.28*5*t))+8).*sin((6.28*25*t)+3.14);
%x = 4*randn(1,M);     % innovation
d = s+x;
% prediction order (FIR filter order + 1)
p = 20
b = zeros(p,1);
% FIR filter order
N  = length(b)-1;
Lx = length(x);
x_init = x(1:500);
Lx_init = length(x_init);
Rx = xcorr(x_init,'unbiased');
disp(Rx(Lx_init));

%alpha = 0.000003;
%alpha = 1/(Rx(Lx_init)*(p+1)*10);
%disp(alpha);

y = zeros(1,Lx);
e = zeros(1,Lx+1);
xx = [zeros(p,1)' x];
rx0 = Rx(Lx_init);

for k = 1:Lx
    if k > Lx_init
      %rx1 = rx0 + ((1/Lx_init)*((xx(k)*xx(k))-(xx(k-Lx_init)*xx(k-Lx_init))));
      %rx1 = rx0 + ((1/Lx_init)*((xx(k)*xx(k))-rx0));
      Rx = xcorr(x(k-Lx_init+1:k),'unbiased');
      rx1 = Rx(Lx_init);
      mu = 1/(rx1*(p+1)*10);
      xn  = xx(k:1:k+N)';
      b = b + 2*mu*e(k)*xn;
      y(k) = xx(k+1:1:k+1+N)*b;
      e(1+k) = d(k) - y(k);
      rx0 = rx1;
    endif
end

e = e(1:Lx);

f = [0:floor((M-1)/2)] / (M*(1/1000));
Y = fft(hamming(M)'.*s');
Y = Y/M;
Y = [Y(1) 2*Y(2:floor((M-1)/2)+1)];

figure(1)
clf
subplot(3,1,1)
plot(t,s,'k')
ylim([-25 25])
xlim([0 (M-1)*T])
title('EKG-Signal','FontName','Arial','FontSize',12)
grid on
subplot(3,1,2)
plot(t,x,'k')
ylim([-25 25])
xlim([0 0.5])
title('50Hz Netzbrummen','FontName','Arial','FontSize',12)
grid on
subplot(3,1,3)
plot(t,d,'k')
ylim([-45 45])
xlim([0 (M-1)*T])
title('EKG-Signal mit �berlagertem Netzbrummen','FontName','Arial','FontSize',12)
grid on

figure(2)
clf
subplot(3,1,1)
plot(t,y,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('Gesch�tztes 50Hz Netzbrummen','FontName','Arial','FontSize',12)
grid on
subplot(3,1,2)
plot(t,d,'k')
ylim([-45 45])
xlim([Lx_init*T (M-1)*T])
title('EKG-Signal mit �berlagertem Netzbrummen','FontName','Arial','FontSize',12)
grid on
subplot(3,1,3)
plot(t,e,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('EKG-Signal st�rbefreit','FontName','Arial', 'FontSize',12)
grid on

figure(3)
clf
subplot(2,1,1)
plot(t,s,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('Original signal','FontName','Arial','FontSize',12)
grid on
subplot(2,1,2)
plot(t,e,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('Noise free signal','FontName','Arial', 'FontSize',12)
grid on

figure(4)
clf
plot(f,abs(Y),'k')
xlim([0 100])
grid on

