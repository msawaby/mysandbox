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
N=2^16;                 % Samples number
Ntransient= 2^12;         % number of simulation points(more accuracy more higher
timex=(0:Ts:(N+Ntransient-1)*Ts)';
%random=[timex 0.257*rand(N+Ntransient,1)-(.257*.5)]; %rpdf
%random=[timex 0*(0.257*rand(N+Ntransient,1)-0.257*rand(N+Ntransient,1))];%tpdf
random=0.5*.225;
Fin_des = 10;
nper = ceil(Fin_des*N/Fs);
% ************************************************************************
% Coefficents of the input sin signal
% ************************************************************************
Fin=nper*Fs/N;			% I+nput signal frequency (Fin = nper*Fs/N)
Ampl=0.6;          	% Input signal amplitude [V]
const=0;
Vin_bias = 0.9;           % Input bias voltage
% About the number of levels in DAC 
levels = 8;
MIN_volt=0;
MAX_volt=1.8;


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
b3=1;
a1=b1;
a2=b2;
c1=0.2817;
c2=5.7104;
r1=b1;
r2=b2;
% echo on;
k=1.38e-23;             % Boltzmann Constant
Temp=300;				% Absolute Temperature in Kelvin
Cf=2.5e-12;				% Integrating Capacitance of the first integrator
alfa=0.999;         	% A=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
                        % A in our exaple =1000 
Amax=1.8;				% Op-amp saturation value [V]
sr=inf;                 % Op-amp slew rate [V/s]
GBW=4000000;			% Op-amp GBW [Hz]
noise1=30e-9;			% 1st int. output noise std. dev. [V/sqrt(Hz)]
b=1;
%echo off;
% Modulator coefficients
Vref=1.8;
finrad=Fin*2*pi;		% Input signal frequency in radians

DAC_step=0.225;
DAC_shift=0.225;

pe = 0.01 / 3;
frac_error = pe * randn(8,1);
disp(frac_error)
UEG0 = DAC_step * ( 1 + frac_error(1) );
UEG1 = DAC_step * ( 1 + frac_error(2) );
UEG2 = DAC_step * ( 1 + frac_error(3) );
UEG3 = DAC_step * ( 1 + frac_error(4) );
UEG4 = DAC_step * ( 1 + frac_error(5) );
UEG5 = DAC_step * ( 1 + frac_error(6) );
UEG6 = DAC_step * ( 1 + frac_error(7) );
UEG7 = DAC_step * ( 1 + frac_error(8) );



% ************************************************************************
% Open Simulink diagram first
% ************************************************************************

options=simset('InitialState', zeros(1,5), 'RelTol', 1e-3, 'MaxStep', 1/Fs );
sim('Sigma_Delta_NLDAC_BiDWA_SL', (N+Ntransient)/Fs, options);	% Starts Simulink simulation

% ************************************************************************
%   Calculates SNR and PSD of the bit-stream and of the signal
% ************************************************************************
w=hann(N)/(N/4);
%echo on;
f=Fin/Fs;			% Normalized signal frequency
fB=N*(bw/Fs);		% Base-band frequency bins

yy12=yout(2+Ntransient:1+N+Ntransient)';
yy_dac = DAC_Output(2+Ntransient:1+N+Ntransient)';
yy_signal = (N/sum(w))*sinusx(yy12(1:N).*w',f,N);
DAC_Distortion = yy_dac ./ DAC_step - yy12 - 1 + yy_signal;
out_bias = sum(yy12)/length(yy12) ;%- sum(DAC_Distortion)/length(DAC_Distortion);
yy1 = yy12-out_bias;
dac_bias = sum(DAC_Distortion) ./ length(DAC_Distortion);
DAC_Distortion = DAC_Distortion - dac_bias;
%echo off;


[snr,ptot,psig,pnoise]=calcSNR(yy1,f,3,fB,w',N);
[snr_dac,ptot_dac,psig_dac,pnoise_dac]=calcSNR(DAC_Distortion,f,1,fB,w',N);

Rbit =(snr-1.76)/6.02;	% Equivalent resolution in bits
Rbit_dac = (snr_dac - 1.76) / 6.02;

% ************************************************************************
% Output Graphs
% ************************************************************************

figure;
clf;
semilogx(linspace(0,Fs/2,N/2), pnoise(1:N/2), 'r');
grid on;
title('PSD of a 2nd-Order Sigma-Delta Modulator')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);


figure;
clf;
semilogx(linspace(0,Fs/2,N/2), psig(1:N/2), 'r');
grid on;
title('PSD of a 2nd-Order Sigma-Delta Modulator')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);

figure;
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
disp(s1)

figure;
clf;
semilogx(linspace(0,Fs/2,N/2), ptot_dac(1:N/2), 'r');
grid on;
title('DAC Distortion')
xlabel('Frequency [Hz]')
ylabel('PSD [dB]')
axis([0 Fs/2 -200 0]);
text(floor(Fs/R),-40, sprintf('SNR = %4.1fdB @ OSR=%d\n',snr_dac,R));
text(floor(Fs/R),-60, sprintf('Rbit = %2.2f bits @ OSR=%d\n',Rbit_dac,R));

s2=sprintf('   SNR(dB)=%1.3f',snr_dac);
s3=sprintf('   Simulation time =%1.3f min',etime(clock,t0)/60);
disp(s2)
disp(s3)