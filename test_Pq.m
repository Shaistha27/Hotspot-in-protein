clc;
close all;
clear all;
% Define the EIIP values for each amino acid
% Define simple discrete sequences
sequence1 = [1, 2, 3, 4];  % Example sequence 1
sequence2 = [4, 3, 2, 1];  % Example sequence 2

% Calculate the Fast Fourier Transform (FFT) of both sequences
RFT1 = fft(sequence1);  % FFT of sequence 1
RFT2 = fft(sequence2);  % FFT of sequence 2

% Calculate the Cross-Spectral Function (Pq)
Pq = abs(RFT1 .* conj(RFT2));  % Magnitude of product of RFT1 and conjugate of RFT2

% Display the result
disp('Cross-Spectral Function (Pq):');
disp(Pq);

