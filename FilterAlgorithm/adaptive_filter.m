clear all;
pkg load signal

% The sampling frequency in Hz.
Fs = 1000;
% Sampling time
T = 1/Fs;
% Buffer length
M = 4000;
% Create a time base
t = (0:M-1)*T;
% Loudspeaker signal
x = (5*sin((6.28*5*t))+5).*(2*sin(6.28*25*t));
% Loudspeaker signal coupled into Microphone (delayed and damped)
x_d = (5*sin(5*(6.28*t-0.3))+5).*(1.5*sin(25*(6.28*t-0.3)));
% Desired signal
d = x_d;
% Predictor order (FIR filter order + 1)
p = 32
% Initialize Coefficients
b = zeros(p,1);
N  = length(b)-1;
% Signal length
Lx = length(x);
% Init length
x_init = x(1:500);
Lx_init = length(x_init);
% Initial cross correlation of Loudspeaker signal
Rx = xcorr(x_init,'unbiased');
rx0 = Rx(Lx_init);
disp(Rx(Lx_init));
% Initialize output and error
y = zeros(1,Lx);
e = zeros(1,Lx+1);

for k = 1:Lx
    % Update coefficients after Initialization
    if k > Lx_init
      % Update correlation value
      rx1 = rx0 + ((1/Lx_init)*((x(k)*x(k))-rx0));
      % Calculate Step size with a factor of 0.5
      beta = 1/(rx1*(p+1)*2);
      % Take input values based on number of coefficients
      xn  = x(k:-1:k-N)';
      % Update filter coefficients
      b = b + 2*beta*e(k)*xn;
      % Predict the desired signal
      y(k) = x(k:-1:k-N)*b;
      % Calculate the difference between desired and predicted signal
      e(1+k) = d(k) - y(k);
      % Buffer old correlation value
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

figure(1)
clf
subplot(2,1,1)
plot(t,x,'k')
ylim([-25 25])
xlim([0.15 1.15])
title('Disturbance','FontName','Arial','FontSize',12)
%grid on
subplot(2,1,2)
plot(t,d,'k')
ylim([-45 45])
xlim([0.15 1.15])
title('Required Signal disturbed','FontName','Arial','FontSize',12)
%grid on

figure(2)
clf
subplot(3,1,1)
plot(t,y,'k')
ylim([-25 25])
xlim([1.15 2.15])
title('Estimated Disturbance','FontName','Arial','FontSize',12)
%grid on
subplot(3,1,2)
plot(t,d,'k')
ylim([-45 45])
xlim([0.15 1.15])
title('Required Signal disturbed','FontName','Arial','FontSize',12)
%grid on
subplot(3,1,3)
plot(t,e,'k')
ylim([-25 25])
xlim([1.15 2.15])
title('Disturbance Free Signal','FontName','Arial', 'FontSize',12)
%grid on

figure(3)
clf
subplot(1,1,1)
plot(t,e,'k')
ylim([-25 25])
xlim([0.15 1.15])
title('Disturbance Free Signal','FontName','Arial', 'FontSize',12)
%grid on

figure(4)
clf
%subplot(2,1,1)
plot(f,abs(Yx),'g')
hold on
plot(f,abs(Ye),'b')
xlim([0 50])
legend('Frequency Response Noise Signal', 'Frequency Response Noise free Signal')
grid on
%{
title('Frequency Response of Disturbance','FontName','Arial', 'FontSize',12)
grid on
subplot(2,1,2)
plot(f,abs(Ye),'g')
ylim([0 0.2])
xlim([0 50])
title('Frequency Response of Disturbance Free Signal','FontName','Arial', 'FontSize',12)
grid on
%}

figure(5)
clf
subplot(2,1,1)
plot(t,x,'k')
ylim([-25 25])
xlim([0.15 1.15])
title('Disturbance','FontName','Arial','FontSize',12)
%grid on
subplot(2,1,2)
plot(t,x_d,'k')
ylim([-25 25])
xlim([0.15 1.15])
title('Disturbance echo','FontName','Arial','FontSize',12)
%grid on
