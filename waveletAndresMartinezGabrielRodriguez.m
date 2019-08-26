clc
clear all 
close all
load('chb19_30_short.mat'); 
samples=808000:840000;
fs=256;
time=samples/fs;
figure(1)
plot(time,x/max(abs(x))),title('Sample EEG')
xlabel('Time,s'),ylabel('Normalized amplitude')



scaled=x/max(abs(x));
wlen=256; %	window	size
h=wlen/4; %hop	size
nSamples=length(samples);
nfft=nSamples; %number	of	fft	samples
K=sum(hamming(wlen,'periodic'))/wlen; %create	the	Hamming	window	and	
		%integrate	the	values
% Create STFT
[s,f,t] = spectrogram(scaled,hamming(wlen,'periodic'),h,nfft,fs,'yaxis');
s = abs(s)/wlen/K;
% Correct	the	DC	and	Nyquist	component
if rem(nfft, 2) % odd nfft excludes Nyquist point
 s(2:end, :) = s(2:end, :).*2;
else % even nfft includes Nyquist point
 s(2:end-1, :) = s(2:end-1, :).*2;
end
% Convert to dB
s = 20*log10(s + 1e-6);
figure(2)
allt = t + time(1);
imagesc(allt, f(1:2200), s(1:2200,:));
axis('xy')
set(gca,'YDir','normal')
set(gca, 'FontName', 'Times New Roman', 'FontSize',20)
xlabel('time, [s]');
ylabel('Frequency, [Hz]');
title('STFT spectrogram of the signal');
colorbar;
h = colorbar;
ylabel(h, 'Magnitude, dB');
colormap jet

%%% CWT
[wt,f]=cwt(scaled,fs/2*pi);
wt=abs(wt)*(2*pi)/fs;
wt = 20*log10(wt + 1e-6);

figure(3)
imagesc(allt, f, wt);
axis('xy')
set(gca,'YDir','normal')
set(gca, 'FontName', 'Times New Roman', 'FontSize',20)
xlabel('time, [s]');
ylabel('Frequency, [Hz]');
title('Amplitude spectrogram of the signal');
colorbar;
h = colorbar;
ylabel(h, 'Magnitude, dB');
colormap jet