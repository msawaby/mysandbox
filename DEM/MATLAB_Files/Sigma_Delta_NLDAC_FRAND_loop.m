clear
clear all
close all
%clc
t0=clock;
% ************************************************************************
% Global varibles
% ************************************************************************
bw=200;  				% Base-band
R=128;                  % Oversampling ratio
Fs=R*2*bw;				% Oversampling frequency
Ts=1/Fs;
N=2^16;                 % Samples number
Ntransient=2^12;         % number of simulation points(more accuracy more higher
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
Ampl=0.380/2;          	% Input signal amplitude [V]
const=0;
Vin_bias=0.9;           % Input bias voltage
% About the number of levels in DAC
levels = 8;
MIN_volt=0;
MAX_volt=1.8;
DAC_step=0.225;

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

avg_loops = 10;
pe_points = 21;
pe = zeros(pe_points,1);
frac_error = zeros(pe_points,avg_loops,8);
snr = zeros(pe_points,1);
snr_dac = zeros(pe_points,1);
snr_tmp = zeros(avg_loops,1);
snr_dac_tmp = zeros(avg_loops,1);
Rbit = zeros(pe_points,1);
Rbit_dac = zeros(pe_points,1);

for i = 1:pe_points,
    pe(i) = 0.0005 * ( i - 1 );
    s1=sprintf('   Element Fraction Error %1.6f.',pe(i));
    disp(s1)
    for j = 1:avg_loops,
        
        s1=sprintf('   Iteration No. %1.0f / %1.0f : %1.0f / %1.0f.',i,pe_points,j,avg_loops);
        disp(s1)
        frac_error(i,j,:) = pe(i) * randn(8,1) / 3;
        %disp(frac_error(i,j,:))
        UEG0 = DAC_step * ( 1 + frac_error(i,j,1) );
        UEG1 = DAC_step * ( 1 + frac_error(i,j,2) );
        UEG2 = DAC_step * ( 1 + frac_error(i,j,3) );
        UEG3 = DAC_step * ( 1 + frac_error(i,j,4) );
        UEG4 = DAC_step * ( 1 + frac_error(i,j,5) );
        UEG5 = DAC_step * ( 1 + frac_error(i,j,6) );
        UEG6 = DAC_step * ( 1 + frac_error(i,j,7) );
        UEG7 = DAC_step * ( 1 + frac_error(i,j,8) );
        
        DAC_shift=0.225;
        
        % ************************************************************************
        % Open Simulink diagram first
        % ************************************************************************
        
        options=simset('InitialState', zeros(1,2), 'RelTol', 1e-3, 'MaxStep', 1/Fs );
        sim('Sigma_Delta_NLDAC_FRAND_SL', (N+Ntransient)/Fs, options);	% Starts Simulink simulation
        
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
        
        
        snr_tmp(j)=calcSNR(yy1,f,3,fB,w',N);
        snr_dac_tmp(j)=calcSNR(DAC_Distortion,f,3,fB,w',N);
        s1=sprintf('   SNR For This Iteration (dB)=%1.3f',snr_tmp(j));
        disp(s1)
        s1=sprintf('   SNR Of DAC(dB)=%1.3f',snr_dac_tmp(j));
        disp(s1)
        
    end
    snr(i) = snr_avg(snr_tmp);
    snr_dac(i) = snr_avg(snr_dac_tmp);
    Rbit(i) =(snr(i)-1.76)/6.02;	% Equivalent resolution in bits
    Rbit_dac(i) = (snr_dac(i) - 1.76) / 6.02;
        
    s1=sprintf('   SNR(dB)=%1.3f',snr(i));
    disp(s1)
    
end


% ************************************************************************
% Output Graphs
% ************************************************************************

figure(1);
clf;
plot(pe, snr, 'r');
grid on;
title('SNR vs. Mismatch')
xlabel('Percentage Mismatch')
ylabel('SNR')

figure(2);
clf;
plot(pe, snr_dac, 'r');
grid on;
title('DAC SNR vs. Mismatch')
xlabel('Percentage Mismatch')
ylabel('SNR')

s2=sprintf('   Simulation time =%1.3f min',etime(clock,t0)/60);
disp(s2)