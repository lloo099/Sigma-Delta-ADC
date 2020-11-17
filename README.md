# Sigma-Delta-ADC
Demonstration of the simulateDSM function, as in the MATLAB Delta Sigma Toolbox, employing its Python port deltasigma.

In the first section, the Noise Transfer Function (NTF) is synthesized for a 2th-order, low-pass modulator.

In each case, the modulator is simulated - with simulateDSM - and its output plotted in terms of time response and spectrum.
The Signal to Noise Ratio (SNR) is evaluated from the spectrum through calculateSNR.
The SNR for different amplitudes is predicted through predictSNR, simulated with simulateSNR and the peak values are extracted with peakSNR.
reference web: http://www.python-deltasigma.io/

The verilog-A models also including in here, you can send e-mail to me if you have any suggestion about it.Thx..

Thomas J.J. Chou address: zhoutomas177@gmail.com
