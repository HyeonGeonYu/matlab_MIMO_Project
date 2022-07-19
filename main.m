%{
  * Original Code
  * https://github.com/Yunseong-Cho/LearningML/blob/master/MLdecoding/Main_ML.m
%}

clear;
clc;
close;
%% Start Main
tic;

example_val ="Ex_5"

%% Ex_1 Maximum Information Rate Design
if example_val == "Ex_1"
  % system parameters
  p_0 = 1;
  Lambda_matrix = diag([300,100,60,30,20]);
  W_matrix = Lambda_matrix;
  B = 5
  [mu,diag_f_matrix] = mu_Cal(Lambda_matrix,W_matrix,p_0,B);
  Gamma_matrix = diag_f_matrix^2*Lambda_matrix;
  
  t_i = transpose([1,2,3,4,5]);
  t_lambda_i = diag(Lambda_matrix);
  t_phi_square_f_i = diag(diag_f_matrix^2);
  t_gamma_i = diag(10*log10(Gamma_matrix)); %dB
  t_rate_i = log2(1+diag(Gamma_matrix));
  t_M_i = 2.^floor(t_rate_i);
  table_result = table(t_i,t_lambda_i,t_phi_square_f_i,t_gamma_i,t_rate_i,t_M_i)

%% Ex_2 Qos Based Design
elseif example_val == "Ex_2"
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
  SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)

  for ss = 1: length(sigma_val)

    R_nn = (sigma_val(ss)^2)*eye(SP.Nr); % Number of receiver antenna
    for expr = 1 : experiments
      
      [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
      [V_matrix,Lambda_matrix] = eig(H'*R_nn^(-1)*H);
      [Lambda_vector,idx_Lm] = sort(diag(Lambda_matrix),'descend');
      Lambda_matrix = Lambda_matrix(idx_Lm,idx_Lm);
      Lambda_matrix = Lambda_matrix(1:B,1:B);
      D_matrix = diag([0.75,0.25]);
      
      gamma_scalar = p_0 / trace( D_matrix * Lambda_matrix^(-1) );
      Gamma_matrix = gamma_scalar*D_matrix;
      tmp = diag(real(Gamma_matrix));
      result_value_tmp_video(ss) = result_value_tmp_video(ss) + tmp(1);
      result_value_tmp_audio(ss) = result_value_tmp_audio(ss) + tmp(2);
    end
    fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
  end

  result_value_tmp_video = 10*log10(result_value_tmp_video/experiments);
  result_value_tmp_audio = 10*log10(result_value_tmp_audio/experiments);
  
  plot(SNR_dB,result_value_tmp_video,'-g+')
  hold on;
  plot(SNR_dB,result_value_tmp_audio,'-r*')
  hold off;
  grid on
  grid minor
  legend('Video Stream','Audio Stream')
  xlabel('Total transmit power/receive  noise(in dB)') 
  ylabel('Sub-channel SNR(in dB)') 

%% Ex_3 (Unweigthed) MMSE Design
elseif example_val == "Ex_3"
  % system parameters
  experiments = 500000
  p_0 = 1;
  B = 2;
  Nt_list = 2:2:6;
  
  SP.Nr = 2;     % Number of receiver antenna
  SNR_dB = 0:2:18;   %(10*log10(p_0/sigma^2))
  sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
  
  result_BER = zeros(length(SNR_dB),length(Nt_list));
  SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)
  W_matrix = eye(B);
  beta_M_i = 3/(4-1); %4QAM
  N_e_i = 2; %4QAM
  
  for Nt_idx = 1: length(Nt_list)
    SP.Nt = Nt_list(Nt_idx);      
    for ss = 1: length(sigma_val)

      R_nn = (sigma_val(ss)^2)*eye(SP.Nr); % Number of receiver antenna
      for expr = 1 : experiments
        [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
        [V_matrix,Lambda_matrix] = eig(H'*R_nn^(-1)*H);
        [Lambda_vector,idx_Lm] = sort(diag(Lambda_matrix),'descend');
        Lambda_matrix = Lambda_matrix(idx_Lm,idx_Lm);
        Lambda_matrix = Lambda_matrix(1:B,1:B);
        V_matrix = V_matrix(:,idx_Lm);
        V_matrix = V_matrix(:,1:B);
        [mu,diag_f_matrix] = mu_Cal(Lambda_matrix,W_matrix,p_0,B);
        Gamma_matrix = W_matrix^(0.5) * mu^(-0.5)*Lambda_matrix^(0.5)-eye(B);
        Gamma_matrix(Gamma_matrix<0)=0;

        tmp = diag(real(Gamma_matrix));
        result_BER(ss,Nt_idx) = result_BER(ss,Nt_idx) + mean(N_e_i*qfunc(sqrt(beta_M_i*tmp)));
      end
      fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
    end
  end

  result_BER = result_BER/experiments;
  
  plot(SNR_dB,result_BER(:,1),'-r')
  hold on
  plot(SNR_dB,result_BER(:,2),'-g+')
  plot(SNR_dB,result_BER(:,3),'-bo')
  hold off
  grid on
  grid minor
  legend('MT = 2','MT = 4','MT = 6')
  xlabel('Total transmit power/receive  noise(in dB)') 
  ylabel('Average BER') 
  ylim([10^(-4) 1]) 
  set(gca, 'YScale', 'log')


%% Ex_4 Equal Error Design
elseif example_val == "Ex_4"
  % system parameters
  experiments = 50000
  p_0 = 1;  
  B_list = 3:1:5;

  SP.Nt = 5;
  SP.Nr = 5;     % Number of receiver antenna
  SNR_dB = 0:2:18;   %(10*log10(p_0/sigma^2))
  sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
  
  result_BER = zeros(length(SNR_dB),length(B_list));
  SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)
  
  beta_M_i = 3/(4-1); %4QAM
  N_e_i = 2; %4QAM
  
  for B_idx = 1: length(B_list)
    B = B_list(B_idx);      
    D_matrix = 1/B*eye(B);
    for ss = 1: length(sigma_val)
      R_nn = (sigma_val(ss)^2)*eye(SP.Nr); % Number of receiver antenna
      for expr = 1 : experiments
        [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
        [V_matrix,Lambda_matrix] = eig(H'*R_nn^(-1)*H);
        [Lambda_vector,idx_Lm] = sort(diag(Lambda_matrix),'descend');
        Lambda_matrix = Lambda_matrix(idx_Lm,idx_Lm);
        Lambda_matrix = Lambda_matrix(1:B,1:B);

        gamma_scalar = p_0/trace(D_matrix*Lambda_matrix^(-1));
        Gamma_matrix = gamma_scalar*D_matrix;
        tmp = real(Gamma_matrix(1));
        result_BER(ss,B_idx) = result_BER(ss,B_idx) + N_e_i*qfunc(sqrt(beta_M_i*tmp));
      end
      fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
    end
  end

  result_BER = result_BER/experiments;
  
  plot(SNR_dB,result_BER(:,1),'-r')
  hold on
  plot(SNR_dB,result_BER(:,2),'-g+')
  plot(SNR_dB,result_BER(:,3),'-bo')
  hold off
  grid on
  grid minor
  legend('B = 3','B = 4','B = 5')
  xlabel('Total transmit power/receive  noise(in dB)') 
  ylabel('Average BER') 
  ylim([10^(-4) 1]) 
  set(gca, 'YScale', 'log')

%% Ex_5 Comparision of Equal-Error and MMSE Design
elseif example_val == "Ex_5"
  % system parameters
  experiments = 50000
  p_0 = 1;
  B = 3;
  SP.Nt = 4;
  SP.Nr = 4;
  SNR_dB = 0:2:18;   %(10*log10(p_0/sigma^2))
  sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
  result_BER = zeros(length(SNR_dB),B+1);
  SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)
  beta_M_i = 3/(4-1); %4QAM
  N_e_i = 2; %4QAM
  
  W_matrix = eye(B);   
  D_matrix = 1/B*eye(B);
  for ss = 1: length(sigma_val)
    R_nn = (sigma_val(ss)^2)*eye(SP.Nr); % Number of receiver antenna
    for expr = 1 : experiments
      [H] = Channel_Gen(SP); % H (channel matrix Nr x Nt)
      [V_matrix,Lambda_matrix] = eig(H'*R_nn^(-1)*H);
      [Lambda_vector,idx_Lm] = sort(diag(Lambda_matrix),'descend');
      Lambda_matrix = Lambda_matrix(idx_Lm,idx_Lm);
      Lambda_matrix = Lambda_matrix(1:B,1:B);
      % MMSE
      [mu,diag_f_matrix] = mu_Cal(Lambda_matrix,W_matrix,p_0,B);
      Gamma_matrix1 = W_matrix^(0.5) * mu^(-0.5)*Lambda_matrix^(0.5)-eye(B);
      Gamma_matrix1(Gamma_matrix1<0)=0;
      tmp1 = transpose(diag(real(Gamma_matrix1)));
      result_BER(ss,1:B) = result_BER(ss,1:B) + N_e_i*qfunc(sqrt(beta_M_i*tmp1));

      % Equal-Error
      gamma_scalar = p_0/trace(D_matrix*Lambda_matrix^(-1));
      Gamma_matrix2 = gamma_scalar*D_matrix;
      tmp2 = real(Gamma_matrix2(1));
      result_BER(ss,B+1) = result_BER(ss,B+1) + N_e_i*qfunc(sqrt(beta_M_i*tmp2));
      
    end
    fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
  end

  result_BER = result_BER/experiments;
  plot(SNR_dB,result_BER(:,1),'-r')
  hold on
  plot(SNR_dB,result_BER(:,2),'-g+')
  plot(SNR_dB,result_BER(:,3),'-bo')
  plot(SNR_dB,result_BER(:,4),'-c')
  hold off
  grid on
  grid minor
  legend('mmse-1','mmse-2','mmse-3','equal-error')
  xlabel('Total transmit power/receive  noise(in dB)') 
  ylabel('Average BER') 
  ylim([10^(-3.8) 1]) 
  set(gca, 'YScale', 'log')
  %}
end