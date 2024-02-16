load('C:\DEM\MATLAB_Files\Loop_NLDAC_DWA.mat','snr','snr_dac','pe');
dwa_snr = snr;
dwa_snr_dac = snr_dac;

load('C:\DEM\MATLAB_Files\Loop_NLDAC_BiDWA.mat','snr','snr_dac','pe');
bidwa_snr = snr;
bidwa_snr_dac = snr_dac;

load('C:\DEM\MATLAB_Files\Loop_NLDAC_FRAND.mat','snr','snr_dac');
frand_snr = snr;
frand_snr_dac = snr_dac;

load('C:\DEM\MATLAB_Files\Loop_NLDAC.mat','snr','snr_dac');
nl_snr = snr;
nl_snr_dac = snr_dac;

pe = pe .* 100;

figure;
plot(pe,dwa_snr,pe,bidwa_snr,pe,nl_snr,pe,frand_snr);
legend('With Data-Weighted Averaging','With Bidirectional Data-Weigthed Aberaging', ...
    'Normal DAC', 'With Full Randamization');
title('SNR vs. Mismatch');
xlabel('(%) Percentage Mismatch');
ylabel('SNR(dB)');

figure;
plot(pe,dwa_snr_dac,pe,bidwa_snr_dac,pe,nl_snr_dac,pe,frand_snr_dac);
legend('With Data-Weighted Averaging','With Bidirectional Data-Weigthed Aberaging', ...
    'Normal DAC', 'With Full Randamization');
title('DAC SNR vs. Mismatch');
xlabel('(%) Percentage Mismatch');
ylabel('SNR(dB)');



load('C:\DEM\MATLAB_Files\Loop_NLDAC_BiDWA_10pc.mat','snr','snr_dac','pe');
bidwa_snr = snr;
bidwa_snr_dac = snr_dac;
pe = pe .* 100;
figure;
plot(pe,bidwa_snr,pe,bidwa_snr_dac);
legend('Modulator SNR','DAC SNR');
title('SNR vs. Mismatch');
xlabel('(%) Percentage Mismatch');
ylabel('SNR(dB)');
