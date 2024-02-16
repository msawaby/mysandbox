function snr = snr_avg(snrs)


snr_linear = 10 .^ ( snrs ./ 10);
avg_snr_linear = sum(snr_linear) ./ length(snr_linear);
snr = dbp(avg_snr_linear);