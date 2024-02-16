clear
clear all
close all
clc
t0=clock;
% ************************************************************************
% Global varibles
% ************************************************************************
bw=200;  				% Base-band
R=128;                  % Oversampling ratio
Fs=R*2*bw;				% Oversampling frequency
Ts=1/Fs;
N=130*500;                 % Samples number
Ntransient=1000;         % number of simulation points(more accuracy more higher
timex=(0:Ts:(N+Ntransient-1)*Ts)';
%random=[timex 0.257*rand(N+Ntransient,1)-(.257*.5)]; %rpdf
%random=[timex 0*(0.257*rand(N+Ntransient,1)-0.257*rand(N+Ntransient,1))];%tpdf
random=0.5*.225;
Fin_des = 100;
nper = ceil(Fin_des*N/Fs);
% ************************************************************************
% Coefficents of the input sin signal
% ************************************************************************
Fin=nper*Fs/N;			% I+nput signal frequency (Fin = nper*Fs/N)
Ampl=0.300;          	% Input signal amplitude [V]
const=0;
Vin_bias=0.9;           % Input bias voltage
% About the number of levels in DAC 
levels = 8;
MIN_volt=0;
MAX_volt=1.8;
DAC_step=0.225;
DAC_shift=0.225;
% ************************************************************************
% Random numbers added to the coefficents (THIS REPRESENTS THE MISMATCH IN
% THE VALUS OF THE CAPACITORS)
% ************************************************************************

% dev1=randn(1,1);
% dev2=randn(1,1);
% dev3=randn(1,1);
% dev4=randn(1,1);
% dev5=randn(1,1);
% dev6=randn(1,1);
% dev7=randn(1,1);

% ************************************************************************
% Coefficents of the structure
% ************************************************************************

% I used them for teseting the effects
% has to be deleted
b1=0.5019;
b2=0.3042;
b3=0;
a1=b1;
a2=b2;
c1=0.2817;
c2=5.7104;
r1=b1;
r2=b2;
% echo on;
k=1.38e-23;			% Boltzmann Constant
Temp=300;				% Absolute Temperature in Kelvin
Cf=2.5e-12;				% Integrating Capacitance of the first integrator
alfa=0.999;         	% A=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
                        % A in our exaple =1000 
Amax=1.8;				    % Op-amp saturation value [V]
sr=inf;                 % Op-amp slew rate [V/s]
GBW=4000000;				% Op-amp GBW [Hz]
noise1=30e-9;			% 1st int. output noise std. dev. [V/sqrt(Hz)]
b=1;
%echo off;
% Modulator coefficients
Vref=1.8;
finrad=Fin*2*pi;		% Input signal frequency in radians

% ************************************************************************
% Open Simulink diagram first
% ************************************************************************

options=simset('InitialState', zeros(1,2), 'RelTol', 1e-3, 'MaxStep', 1/Fs );
sim('Sigma_Delta_Linear_DAC_SL', (N+Ntransient)/Fs, options);	% Starts Simulink simulation

% ************************************************************************
%   Calculates SNR and PSD of the bit-stream and of the signal
% ************************************************************************
w=hann(N)/(N/4);
%echo on;
f=Fin/Fs;			% Normalized signal frequency
fB=N*(bw/Fs);		% Base-band frequency bins
yy12=yout(2+Ntransient:1+N+Ntransient)';
out_bias = sum(yy12)/length(yy12);
yy1 = yy12-out_bias;
%echo off;

%ptot=zeros(1,N);
% spect = fft(yy1.*w, N);
% loglog(abs(spect));
[snr,ptot,psig,pnoise]=calcSNR(yy1,f,3,fB,w',N);
Rbit=(snr-1.76)/6.02;	% Equivalent resolution in bits

% ************************************************************************
% Output Graphs
% ************************************************************************

figure(1);
clf;
semilogx(linspace(0,Fs/2,N/2), pnoise(1:N/2), 'r');
grid on;
title('PSD of a 2nd-Order Sigma-Delta Modulator')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);


figure(2);
clf;
semilogx(linspace(0,Fs/2,N/2), psig(1:N/2), 'r');
grid on;
title('PSD of a 2nd-Order Sigma-Delta Modulator')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);

figure(3);
clf;
semilogx(linspace(0,Fs/2,N/2), ptot(1:N/2), 'r');
hold on;
title('PSD of a 2nd-Order Sigma-Delta Modulator (detail)')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);
grid on;
hold off;
text(floor(Fs/R),-40, sprintf('SNR = %4.1fdB @ OSR=%d\n',snr,R));
text(floor(Fs/R),-60, sprintf('Rbit = %2.2f bits @ OSR=%d\n',Rbit,R));

s1=sprintf('   SNR(dB)=%1.3f',snr);
s2=sprintf('   Simulation time =%1.3f min',etime(clock,t0)/60);
disp(s1)
disp(s2)