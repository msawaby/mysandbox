function out = slew(in,alfa,sr,GBW,Ts)
% Models the operational amplifier finite bandwidth and slew rate
% for a discrete time integrator (by S. Brigati, P. Malcovati)
%
% out = slew(in,alfa,sr,GBW,Ts)
%
% in:		Input signal amplitude
% alfa:		Effect of finite gain (ideal amplifier alfa=1)
% sr:		Slew rate in V/s
% GBW:		Gain-bandwidth product of the integrator loop gain in Hz
% Ts:		Sample time in s
%
% out:		Output signal amplitude

tau=1/(2*pi*GBW);  % Time constant of the integrator
Tmax = Ts/2;

slope=alfa*abs(in)/tau;

if slope > sr			% Op-amp in slewing
	
	tsl = abs(in)*alfa/sr - tau;  % Slewing time
	
	if tsl >= Tmax
		error = abs(in) - sr*Tmax;
	else
		texp = Tmax - tsl;
		error = abs(in)*(1-alfa) + (alfa*abs(in) - sr*tsl) * exp(-texp/tau);
	end
	
else					% Op-amp in linear region
	texp = Tmax;
	error = abs(in)*(1-alfa) + alfa*abs(in) * exp(-texp/tau);
end

out = in - sign(in)*error;
