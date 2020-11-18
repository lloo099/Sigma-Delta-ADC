function [] = adc_Agilent16902a(sPath,nLen,nStart,fSampling)
% Example program routine to generate FFT plots and determine the dynamic performance of
% a dataconverter from the data records taken with a Cadence simulation and Logic Analyzer
% System(Agilent 16902A). 
% Data was read into the following MATLAB program routine. 

%       adc_Agilent16902a(sPath,nLen,nStart,fSampling)
% Example: adc_Agilent16902a('CHIP0_CH1_327680_106.csv',2097151,1,32768/20)
% Example: adc_Agilent16902a('CHIP0_CH1_19660800_9915.csv',65536,1,19660800/20)
% Example: adc_Agilent16902a('CHIP0_CH1_19660800_99915.csv',65536,1,19660800/20)
% Example: adc_Agilent16902a('CHIP0_CH1_19660800_199965.csv',65536,1,19660800/20)



%Start ADC Dynamic Performance Test Routine
close all;

%%%**********************************************************
%%%******* USER INPUT PARAMETERS ****************************
%%%**********************************************************
%%% 1) Simulation result
%filename = './sim_chipcal32_ww.csv';  %% fiename
%SiliconResult=0;  % Result from silicon = 1 or simulation = 0

%%% 2) Silicon result - 1
filename = sPath;  %% fiename

SiliconResult=1;  % Result from silicon = 1 or simulation = 0

%%% 3) Silicon result - 2
%filename = './19M_109703_3p0.csv';  %% fiename
%SiliconResult=1;  % Result from silicon = 1 or simulation = 0

fclk=fSampling;         % input clock frequency
numbit=14;              % number of point
%**********************************************************

xmax=fclk/2;
xmin=1;
ymax=0;

if (SiliconResult == 1) 
    numpt=nLen      % No of FFT
    fid=fopen(filename,'r');
    fgetl(fid); % remove the first title line
    [v1,count]=fscanf(fid,'%f,%x,%f %*s',[3,numpt]); 
    fclose(fid);
    v1=v1';
    code=v1(:,2);


else 
    numpt=64;      % No of FFT
    fid=fopen(filename,'r');
    fgetl(fid); % remove the first title line
    [v1,count]=fscanf(fid,'%f, %f', [2 numpt]); 
    fclose(fid);     
    v1=v1';
    code=v1(:,2)*2^numbit; %auto scaling
end

if (max(code)==2^numbit-1) | (min(code)==0)
    disp('Warning: ADC may be clipping!!!'); 
end

%Recenter the digital sine wave
Dout=code-(2^numbit-1)/2;
if (SiliconResult == 1) 
    Doutw=Dout.*hanning(numpt);
else
    Doutw=Dout;
end

figure;
plot(1:numpt,Dout);

%Performing the Fast Fourier Transform 
Dout_spect=fft(Doutw); 
Dout_spect(1)=min(Dout_spect(2:numpt/2));

%Recalculate to dB
Dout_dB=20*log10(abs(Dout_spect));

%Display the results in the frequency domain with an FFT plot
figure;
maxdB=max(Dout_dB(2:numpt/2))

plot([0:numpt/2-1].*fclk/(numpt),Dout_dB(1:numpt/2)-maxdB);
grid on; 
title('FFT PLOT'); 
xlabel('ANALOG INPUT FREQUENCY (Hz)');
ylabel('AMPLITUDE (dB)');

%Calculate SNR, SINAD, THD and SFDR values
%Find the signal bin number, DC = bin 1
fin=find(Dout_dB(2:numpt/2)==maxdB) 

%Span of the input frequency on each side 
span=max(round(numpt/200),3) 
%span=5;

%Approximate search span for harmonics on each side 
spanh=2;
%Determine power spectrum
spectP=(abs(Dout_spect)).*(abs(Dout_spect)); 
%Find DC offset power 
Pdc=sum(spectP(1:span)); 
%Extract overall signal power 
Ps=sum(spectP(fin-span:fin+span));
%Vector/matrix to store both frequency and power of signal and harmonics
Fh=[]; 
%The 1st element in the vector/matrix represents the signal, the next element represents
%the 2nd harmonic, etc.
Ph=[]; 

%Find harmonic frequencies and power components in the FFT spectrum 
for har_num=1:10
  %Input tones greater than fSAMPLE are aliased back into the spectrum
  tone=rem((har_num*(fin-1)+1)/numpt,1); 
  if tone>0.5 
     %Input tones greater than 0.5*fSAMPLE (after aliasing) are reflected
     tone=1-tone;
  end 
  Fh=[Fh tone]; 
  %For this procedure to work, ensure the folded back high order harmonics do not overlap 
  %with DC or signal or lower order harmonics 
  har_peak=max(spectP(round(tone*numpt)-spanh:round(tone*numpt)+spanh)); 
  har_bin=find(spectP(round(tone*numpt)-spanh:round(tone*numpt)+spanh)==har_peak);
  har_bin=har_bin+round(tone*numpt)-spanh-1;
  Ph=[Ph sum(spectP(har_bin-1:har_bin+1))];
end

%Determine the total distortion power 
Pd=sum(Ph(2:5)); 
%Determine the noise power 
Pn=sum(spectP(1:numpt/2))-Pdc-Ps-Pd;

%format;
A=(max(code)-min(code))/2^numbit ;
AdB=20*log10(A);
SINAD=10*log10(Ps/(Pn+Pd));
SNR=10*log10(Ps/Pn);
disp('THD is calculated from 2nd through 5th order harmonics');
THD=10*log10(Pd/Ph(1));
SFDR=10*log10(Ph(1)/max(Ph(2:10)));
disp('Signal & Harmonic Power Components:');
HD=10*log10(Ph(1:10)/Ph(1));
ENOB=(SINAD - 1.763)/6.02;


text(0.15*(xmax-xmin)/1, ymax-10, sprintf('fin = %.5g Hz', fin*fclk/numpt),'Fontweight','bold','Fontsize',12);
text(0.45*(xmax-xmin)/1, ymax-10, sprintf('SNDR = %.1f dB', SINAD),'Fontweight','bold','Fontsize',12);
text(0.75*(xmax-xmin)/1, ymax-10, sprintf('ENOB = %.1f', ENOB),'Fontweight','bold','Fontsize',12); 
text(0.15*(xmax-xmin)/1, ymax-30, sprintf('SNR = %.5g dB', SNR),'Fontweight','bold','Fontsize',12);
text(0.45*(xmax-xmin)/1, ymax-30, sprintf('THD = %.1f dB', THD),'Fontweight','bold','Fontsize',12);
text(0.75*(xmax-xmin)/1, ymax-30, sprintf('SFDR = %.1f dB', SFDR),'Fontweight','bold','Fontsize',12); 


% histogram boundaries
minbin=min(Dout);
maxbin=max(Dout);
% histogram
h = hist(Dout, minbin:maxbin);
% cumulative histogram
ch = cumsum(h);

% transition levels
T = -cos(pi*ch/sum(h));
% linearized histogram
hlin = T(2:end) -T(1:end-1);
% truncate at least first and last
trunc=2;
hlin_trunc = hlin(1+trunc:end-trunc);
% calculate lsb size and dnl
lsb= sum(hlin_trunc) / (length(hlin_trunc));
dnl= [0 hlin_trunc/lsb-1];
% calculate inl
inl= cumsum(dnl);

inl_max = max(inl);
inl_min = min(inl);
inl_offset = (inl_max+inl_min)/2
inl = inl-inl_offset;

figure;
subplot(2,1,1);
plot(dnl);
axis tight;
grid on;
title('DNL vs. CODE');
xlabel('DIGITAL OUTPUT CODE');
ylabel('DNL (LSB)');
subplot(2,1,2);
plot(inl);
axis tight;
grid on;
title('INL vs. CODE');
xlabel('DIGITAL OUTPUT CODE');
ylabel('INL(LSB)');


