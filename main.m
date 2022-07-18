%{
  * Original Code
  * https://github.com/Yunseong-Cho/LearningML/blob/master/MLdecoding/Main_ML.m
%}

clear;
clc;
close;

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

  experiments = 10000

  p_0 = 1;
  B = 2;
  SP.Nt = 3;      % Number of transmitter antenna
  SP.Nr = 3;     % Number of receiver antenna
  SNR_dB = linspace(0,20,11);   %(10*log10(p_0/sigma^2))
  sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
  result_value_tmp_video = zeros(size(SNR_dB));
  result_value_tmp_audio = zeros(size(SNR_dB));

  for ss = 1: length(sigma_val)
    for expr = 1 : experiments
      R_nn = sigma_val(ss)*eye(3); % Number of receiver antenna
      SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)
      [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
      [V_matrix,Lambda_matrix] = eig(H'*R_nn^(-1)*H);
      [Lambda_vector,idx_Lm] = sort(diag(Lambda_matrix),'descend');
      Lambda_matrix = Lambda_matrix(idx_Lm,idx_Lm);
      Lambda_matrix = Lambda_matrix(1:B,1:B);
      V_matrix = V_matrix(:,idx_Lm);
      V_matrix = V_matrix(:,1:B);
      D_matrix = diag([0.76,0.24]);
      gamma_scalar = p_0 / trace( D_matrix * Lambda_matrix^(-1) );
      Gamma_matrix = gamma_scalar*D_matrix;
      tmp = diag(real(Gamma_matrix));
    end
    result_value_tmp_video(ss) = result_value_tmp_video(ss) + tmp(1);
    result_value_tmp_audio(ss) = result_value_tmp_audio(ss) + tmp(2);
  end
  result_value_tmp_video = result_value_tmp_video/experiments;
  result_value_tmp_audio = result_value_tmp_audio/experiments;
  plot(SNR_dB,result_value_tmp_video)
  hold on;
  plot(SNR_dB,result_value_tmp_audio)
  hold off;
end