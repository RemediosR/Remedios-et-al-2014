% fft_ex.m
%
% This is an example of computing the FFT of 
% the addition of two cos functions. 
%
% Written by Raymond Roberts and Graham Roberts,
% students at Curtin University of Technology
% September 1996

clc

fprintf('**************************************************************\n\n');
fprintf('This file illustrates how to perform the FFT on a time\n');
fprintf('signal using MATLAB\n\n\n');
fprintf('In this case, the time signal is the sum of two cosines\n\n');

echo on
% Sample at 10 Hz, and use a 500 second long sample -> 5001 samples
t = 0:.1:500;

% Calculate the sum of the two cosines
x = cos(3*2*pi*t) + cos(1*2*pi*t);

% Take an 8192-point FFT (power of 2 greater than 5000)
f = fft(x,8192);

% Construct a frequency axis
Freq = -5:10/8192:5-1/8192;

% Plot frequency, magnitude.  fftshift centers around zero.
plot(Freq, abs(fftshift(f)));

echo off