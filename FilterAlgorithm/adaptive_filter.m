clear all;
pkg load signal

% Sampling rate
Fs = 1000;

T = 1/Fs;
M = 4000;
t = (0:M-1)*T;

s = 4*randn(1,M);
x = (5*sin((6.28*5*t))+5).*(2*sin(6.28*25*t));
x_d = (5*sin(5*(6.28*t-0.3))+5).*(1.5*sin(25*(6.28*t-0.3)));

d = s+x_d;
% prediction order (FIR filter order + 1)
p = 32
b = zeros(p,1);
% FIR filter order
N  = length(b)-1;
Lx = length(x);
x_init = x(1:500);
Lx_init = length(x_init);
Rx = xcorr(x_init,'unbiased');
disp(Rx(Lx_init));

y = zeros(1,Lx);
e = zeros(1,Lx+1);
rx0 = Rx(Lx_init);

for k = 1:Lx
    if k > Lx_init
      rx1 = rx0 + ((1/Lx_init)*((x(k)*x(k))-rx0));
      beta = 1/(rx1*(p+1)*10);
      xn  = x(k:-1:k-N)';
      b = b + 2*beta*e(k)*xn;
      y(k) = x(k:-1:k-N)*b;
      e(1+k) = d(k) - y(k);
      rx0 = rx1;
    endif
end

e = e(1:Lx);

f = [0:floor((M-1)/2)] / (M*(1/1000));

Yx = fft(hamming(M)'.*x_d');
Yx = Yx/M;
Yx = [Yx(1) 2*Yx(2:floor((M-1)/2)+1)];

Ye = fft(hamming(M)'.*e');
Ye = Ye/M;
Ye = [Ye(1) 2*Ye(2:floor((M-1)/2)+1)];

Ys = fft(hamming(M)'.*s');
Ys = Ys/M;
Ys = [Ys(1) 2*Ys(2:floor((M-1)/2)+1)];

figure(1)
clf
subplot(3,1,1)
plot(t,s,'k')
ylim([-25 25])
xlim([0 (M-1)*T])
title('Required Signal','FontName','Arial','FontSize',12)
grid on
subplot(3,1,2)
plot(t,x,'k')
ylim([-25 25])
xlim([0 0.5])
title('Disturbance','FontName','Arial','FontSize',12)
grid on
subplot(3,1,3)
plot(t,d,'k')
ylim([-45 45])
xlim([0 (M-1)*T])
title('Required Signal disturbed','FontName','Arial','FontSize',12)
grid on

figure(2)
clf
subplot(3,1,1)
plot(t,y,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('Estimated Disturbance','FontName','Arial','FontSize',12)
grid on
subplot(3,1,2)
plot(t,d,'k')
ylim([-45 45])
xlim([Lx_init*T (M-1)*T])
title('Required Signal disturbed','FontName','Arial','FontSize',12)
grid on
subplot(3,1,3)
plot(t,e,'k')
ylim([-25 25])
xlim([Lx_init*T (M-1)*T])
title('Disturbance Free Signal','FontName','Arial', 'FontSize',12)
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
title('Disturbance Free Signal','FontName','Arial', 'FontSize',12)
grid on

figure(4)
clf
subplot(3,1,1)
plot(f,abs(Yx),'k')
xlim([0 200])
title('Frequency Response of Disturbance','FontName','Arial', 'FontSize',12)
subplot(3,1,2)
plot(f,abs(Ye),'k')
ylim([0 0.2])
xlim([0 200])
title('Frequency Response of Disturbance Free Signal','FontName','Arial', 'FontSize',12)
subplot(3,1,3)
plot(f,abs(Ys),'k')
ylim([0 0.2])
xlim([0 200])
title('Frequency Response of Required Signal','FontName','Arial', 'FontSize',12)
grid on

figure(5)
clf
subplot(2,1,1)
plot(t,x,'k')
ylim([-25 25])
xlim([0 0.5])
title('Disturbance','FontName','Arial','FontSize',12)
grid on
subplot(2,1,2)
plot(t,x_d,'k')
ylim([-25 25])
xlim([0 0.5])
title('Disturbance echo','FontName','Arial','FontSize',12)
grid on
