
clear;
clc;
close;
%% Start Main
tic;

example_val ="Ex_1"

%% Ex_1 Maximum Information Rate Design
if example_val == "Ex_1"
  % system parameters
  experiments = 100000;
  E_tr = 2;

  B = 2;  
  SP.Nt = 2;      % Number of transmitter antenna
  SP.Nr = 2;     % Number of receiver antenna
  SNR_dB = -10:5:30;   %(10*log10( (E_tr/B) / (tr(Rn_matrix)/M) ))
  sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
  result_mse = zeros(6,length(SNR_dB));
  SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)

  for ss = 1: length(sigma_val)


    for expr = 1 : experiments
      
      [H] = Channel_Gen(SP); 
      
      J_rx_matrix = 1/(sigma_val(ss)^2) * H' * H;
      J_tx_matrix = 1/(sigma_val(ss)^2) * H * H';
      
      result_mse_RxMF = 2 - (trace(J_rx_matrix))^2 / (trace(J_rx_matrix^2+J_rx_matrix));
      result_mse_RxZF = 2 - 4 /( 2+ trace(J_rx_matrix^-1));
      result_mse_RxWF = 2 - trace((J_rx_matrix+eye(2))^-1*J_rx_matrix);

      result_mse_TxMF = 2 - (trace(J_tx_matrix))^2 / (trace(J_tx_matrix^2+J_tx_matrix));
      result_mse_TxZF = 2 - 4 /( 2+ trace(J_tx_matrix^-1));
      result_mse_TxWF = 2 - trace((J_tx_matrix+eye(2))^-1*J_tx_matrix);

      result_mse(1,ss) = result_mse(1,ss) + result_mse_RxMF;
      result_mse(2,ss) = result_mse(2,ss) + result_mse_RxZF;
      result_mse(3,ss) = result_mse(3,ss) + result_mse_RxWF;

      result_mse(4,ss) = result_mse(4,ss) + result_mse_TxMF;
      result_mse(5,ss) = result_mse(5,ss) + result_mse_TxZF;
      result_mse(6,ss) = result_mse(6,ss) + result_mse_TxWF;
    end
    fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
  end
  result_mse = real(result_mse/experiments);
  
  plot(SNR_dB,result_mse(1,:),'-r|')
  hold on
  plot(SNR_dB,result_mse(2,:),'-go')
  plot(SNR_dB,result_mse(3,:),'-b*')
  plot(SNR_dB,result_mse(4,:),'-c^')
  plot(SNR_dB,result_mse(5,:),'-mv')
  plot(SNR_dB,result_mse(6,:),'-yd')
  hold off
  grid on
  grid minor
  legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')
  xlabel('Es/N0 in dB') 
  ylabel('MSE') 
  ylim([10^(-2) 2]) 
  set(gca, 'YScale', 'log')

elseif example_val == "Ex_2"
    % system parameters
    experiments = 1000;
    E_tr = 2;
    B = 2;  
    SP.Nt = 2;      % Number of transmitter antenna
    SP.Nr = 2;     % Number of receiver antenna
    SNR_dB = -10:5:30;   %(10*log10( (E_tr/B) / (tr(Rn_matrix)/M) ))
    sigma_val = (10 .^ (SNR_dB / 10) ) .^ (-0.5);
    result_mse = zeros(6,length(SNR_dB));
    SP.H_type = 'Rayleigh'; % Channel type (Rayleigh or ...)
  
    for ss = 1: length(sigma_val)
  
  
      for expr = 1 : experiments
        
        [H] = Channel_Gen(SP); 
        
        J_rx_matrix = 1/(sigma_val(ss)^2) * H' * H;
        J_tx_matrix = 1/(sigma_val(ss)^2) * H * H';
        
        result_mse_RxMF = 2 - (trace(J_rx_matrix))^2 / (trace(J_rx_matrix^2+J_rx_matrix));
        result_mse_RxZF = 2 - 4 /( 2+ trace(J_rx_matrix^-1));
        result_mse_RxWF = 2 - trace((J_rx_matrix+eye(2))^-1*J_rx_matrix);
  
        result_mse_TxMF = 2 - (trace(J_tx_matrix))^2 / (trace(J_tx_matrix^2+J_tx_matrix));
        result_mse_TxZF = 2 - 4 /( 2+ trace(J_tx_matrix^-1));
        result_mse_TxWF = 2 - trace((J_tx_matrix+eye(2))^-1*J_tx_matrix);
  
        result_mse(1,ss) = result_mse(1,ss) + result_mse_RxMF;
        result_mse(2,ss) = result_mse(2,ss) + result_mse_RxZF;
        result_mse(3,ss) = result_mse(3,ss) + result_mse_RxWF;
  
        result_mse(4,ss) = result_mse(4,ss) + result_mse_TxMF;
        result_mse(5,ss) = result_mse(5,ss) + result_mse_TxZF;
        result_mse(6,ss) = result_mse(6,ss) + result_mse_TxWF;
      end
      fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',ss,toc,(toc/3600))
    end
    result_mse = real(result_mse/experiments);
    
    plot(SNR_dB,result_mse(1,:),'-r|')
    hold on
    plot(SNR_dB,result_mse(2,:),'-go')
    plot(SNR_dB,result_mse(3,:),'-b*')
    plot(SNR_dB,result_mse(4,:),'-c^')
    plot(SNR_dB,result_mse(5,:),'-mv')
    plot(SNR_dB,result_mse(6,:),'-yd')
    hold off
    grid on
    grid minor
    legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')
    xlabel('Es/N0 in dB') 
    ylabel('MSE') 
    ylim([10^(-2) 2]) 
    set(gca, 'YScale', 'log')

end