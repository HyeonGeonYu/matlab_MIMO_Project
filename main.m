%%
%{

  * Original Code
  * https://github.com/Yunseong-Cho/LearningML/blob/master/MLdecoding/Main_ML.m
  
%}

clear;
clc;

% system parameters

SP.Nt = 2;      % Number of transmitter antenna
SP.Nr = 2;     % Number of receiver antenna
SP.SNR_dB = linspace(-8,8,5);   % s_hat = G*H*F*s + G*n  
                                % s_hat (received vector B x 1)
                                % s (transmitted vector B x 1) 
                                % n (noise vector(Nr x 1)
                                % G (decoder matrix B x Nr)
                                % F (precoder matrix Nt x B)
                                % snr = p/N0
                                        
SP.H_type = 'Rayleigh';         % Channel type (Rayleigh or ...)

[H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)


