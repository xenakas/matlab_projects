function clptdscf(f,F,xunit)
% Plot a discrete spectrum [f F]
% CALL: clptdscf(f,F,xunit)
%  Input to function is 
%	f 	- frequencies 
%	F 	- spectral values
%	xunit	- units of frequency (Hz or rad/sec)
% Compute plotting arguments for axis
nl=length(f);		% Number of f points
fmin=min(f);
fmax=max(f);
Fmax=max(F);
% Plotting range, lengthen axes by 10%  
Fmaxp=Fmax+.1*Fmax;
fminp=fmin-.1*fmax;
fmaxp=fmax+.1*fmax;
%
axis([fminp fmaxp 0 Fmaxp])   % Manual scaling
for I=1:nl,
fplots=[f(I) f(I)];
Fplot=[0 F(I)];
plot(fplots,Fplot)	      % Plot one line at a time 	
hold on
end
title('Frequency Spectrum ')
ylabel('Amplitude')
xlabel(['Frequency in ', xunit])

