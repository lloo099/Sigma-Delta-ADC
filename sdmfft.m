order=2; 
OSR=256; 
nLev =3; 
opt=3; 
Xlim=0.9;
f=0;
H = synthesizeNTF(order,OSR,0,opt); 
form = 'CIFB';
[a,g,b,c] = realizeNTF (H,form);
b(2:end) = 0;
ABCD = stuffABCD (a,g,b,c,form); 
[ABCDs,umax]=scaleABCD (ABCD,nLev,f,Xlim)
[a,g,b,c]= mapABCD(ABCDs,form) %scaled ABCD matrix 
[Ha, Ga] =calculateTF (ABCDs)
f =linspace (0,0.5,1000);
z =exp (2i*pi*f);
magHa =dbv(evalTF(Ha,z));
magGa = dbv(evalTF(Ga,z));

Nfft=16384;
fbin=31;
t=(0:Nfft-1);
%u=0.5*sin(2*pi*tone_bin/Nfft*t);
v=csvread('/Users/jiajunzhou/Documents/MATLAB/Examples/R2019a/signal/ChebyshevTypeIFilterDesignExample/sdadc16384_02.csv',1,1);

%spec=fft(v)/(Nfft*(nLev-1)/2); 
spec=fft(v.*hann(Nfft))/(Nfft*(nLev-1)/4);
ntf0 = synthesizeNTF(2,OSR,1);
k = mean(abs(v)/mean(v.^2))
ntf = ntf0 / (k +(1-k)*ntf0);
nb = 3;
% Compute windowed FFT and NBW 
w = hann(Nfft);
w1 = norm(w,1);
w2 = norm(w,2);
NBW = (w2/w1)^2 
V=fft(w.*v)/(w1/2);
% Compute SNR
signal_bins = fbin + [-(nb-1)/2:(nb-1)/2];
inband_bins = 0:Nfft/(2*OSR);
noise_bins = setdiff(inband_bins,signal_bins);
snr = dbp( sum(abs(V(signal_bins+1)).^2) / sum(abs(V(noise_bins+1)).^2))

f=linspace (0,0.5,(Nfft/2)+1); 
Sqq=4*(evalTF(Ha,exp(2i*pi*f))/(nLev-1).^2/3);
subplot(2,1,1);
plot(f, dbv(spec(1 :(Nfft/2)+1)), 'b');
title('Spectrum PLOT of scaled DSM with hann window');
xlabel ('normalized frequency');
ylabel('dBFS');
hold on;
subplot(2,1,2);
plot(f,dbv(Sqq*NBW),'m','Linewidth',1);
ylabel('dBFS');
xlabel ('normalized frequency');
text_handle = text(0.05,-40, sprintf('SNR = %4.1fdB @ OSR= %d',snr,OSR),'vert','middle');
text(0.5, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
