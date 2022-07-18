%{
  * Original Code
  * https://github.com/Yunseong-Cho/LearningML/blob/master/MLdecoding/Main_ML.m
%}

clear;
clc;

example_val ="Ex_2"

%% Ex_1
if example_val == "Ex_1"
  % system parameters
  p_0 = 1;
  Lambda_matrix = diag([300,100,60,30,20]);
  W_matrix = Lambda_matrix;
  mu = mu_Cal(Lambda_matrix,W_matrix,p_0);
  diag_f_matrix = ( eye(5) / (mu^(0.5)) - (Lambda_matrix)^(-1))^(0.5);
  diag_f_matrix(diag_f_matrix<0) = 0;
  Gamma_matrix = diag_f_matrix^2*Lambda_matrix;
  
  t_i = transpose([1,2,3,4,5]);
  t_lambda_i = diag(Lambda_matrix);
  t_phi_square_f_i = diag(diag_f_matrix^2);
  t_gamma_i = diag(10*log10(Gamma_matrix)); %dB
  t_rate_i = log2(1+diag(Gamma_matrix));
  t_M_i = 2.^floor(t_rate_i);
  patients = table(t_i,t_lambda_i,t_phi_square_f_i,t_gamma_i,t_rate_i,t_M_i)

%% Ex_2 
elseif example_val == "Ex_2"
  %% Ex_2
  % system parameters
  SP.Nt = 3;      % Number of transmitter antenna
  SP.Nr = 3;     % Number of receiver antenna
  B = 2;
  SP.SNR_dB = linspace(-8,8,5);   

                                          
  SP.H_type = 'Rayleigh';         % Channel type (Rayleigh or ...)

  [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
  
end