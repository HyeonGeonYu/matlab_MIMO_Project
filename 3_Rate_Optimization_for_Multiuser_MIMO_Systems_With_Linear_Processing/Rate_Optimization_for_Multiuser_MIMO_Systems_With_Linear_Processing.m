clear;
addpath C:\Users\Hyeongeon\ggplab;
global QUIET
QUIET = 1;
User_Num_K = 2;
Tx_antenna_num_N_T = 4;
Rx_antenna_num_N_R_k = 2;
each_user_data_num_M_k = 2 * ones(1, User_Num_K);
%each_user_data_num_M_k = transpose(randi([1, 4], 1, User_Num_K));
total_data_num_N_d = sum(each_user_data_num_M_k);


% 시뮬레이션 세팅
SNR_values = 10:1:10; % SNR 값 범위
num_simulations = 1; % 전체 시뮬레이션 수
num_iteration = 30; % iteration 수
sigma_n = 1;
%P_max = 1;
for snr_index = 1:length(SNR_values)
    % sigma_n = 1/sqrt(10^(SNR/10));
    SNR = SNR_values(snr_index);
    P_max = 10^(SNR/10);
    num_errors = 0;

    for sim = 1:num_simulations
        %transmit_data_symbol_d = sqrt(1/2) * (1 * randn(total_data_num_N_d,1) + 1j * randn(total_data_num_N_d,1));
        %noise_n = sigma_n *sqrt(1/2)* (1 * randn(total_data_num_N_d, 1) + 1j * randn(total_data_num_N_d, 1));       
        each_user_channel_coefficient_H = sqrt(1/2) * (1 * randn(Tx_antenna_num_N_T,User_Num_K*Rx_antenna_num_N_R_k) + 1j * randn(Tx_antenna_num_N_T,User_Num_K*Rx_antenna_num_N_R_k)); 
        
        %% 0) initialization
        user_filter_T = []; % UL-transmitter filter(DL-receiver filter)
        for tmp_k = 1: User_Num_K
            [U,S,V] = svd(each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k));
            user_filter_T = blkdiag(user_filter_T,V');
        end
        % calcuation of the UL-receiver filter(DL-transmitter filter) and beta
        % UL- power allocation
        power_allocation_P = P_max/total_data_num_N_d * eye(total_data_num_N_d);
        % equation (4), UL-receiver filter(DL-transmitter filter)
        bs_filter_U = [];
        scailed_factor_beta = [];
        for tmp_k = 1: User_Num_K
            U_k_beta_k = inv(each_user_channel_coefficient_H * user_filter_T * power_allocation_P * user_filter_T' * each_user_channel_coefficient_H' + sigma_n^2*eye(Tx_antenna_num_N_T)) * each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k)  *  user_filter_T(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k) * power_allocation_P(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k * Rx_antenna_num_N_R_k,1+(tmp_k-1) * (Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k); 
            beta_k = diag(sqrt(sum(abs(U_k_beta_k).^2)));
            U_k = U_k_beta_k*inv(beta_k);
            bs_filter_U = [bs_filter_U,U_k];
            scailed_factor_beta = blkdiag(scailed_factor_beta,beta_k);
        end 

        %% iteration
        for iter_tmp = 1:num_iteration
            %% 1) Optimization in the uplink channel
            % equation (7)
            diagnal_matrix_D = [];
            for tmp_i = 1 : total_data_num_N_d 
                diagnal_matrix_D_ii = scailed_factor_beta(tmp_i,tmp_i)^(2) * ( bs_filter_U(:,tmp_i)' * each_user_channel_coefficient_H * user_filter_T(:,tmp_i) * user_filter_T(:,tmp_i)' * each_user_channel_coefficient_H' * bs_filter_U(:,tmp_i)) - 2 * scailed_factor_beta(tmp_i,tmp_i) * real(bs_filter_U(:,tmp_i)' * each_user_channel_coefficient_H * user_filter_T(:,tmp_i)) + 1;
                diagnal_matrix_D = [diagnal_matrix_D,diagnal_matrix_D_ii]; 
            end
            diagnal_matrix_D = real(diag(diagnal_matrix_D));
            % equation (8)
            psi_matrix = zeros(total_data_num_N_d);
            for tmp_i = 1 : total_data_num_N_d 
                for tmp_j = 1 : total_data_num_N_d 
                    if tmp_i ~= tmp_j
                        psi_matrix(tmp_i,tmp_j) = real(bs_filter_U(:,tmp_i)' * each_user_channel_coefficient_H * user_filter_T(:,tmp_j) * user_filter_T(:,tmp_j)' * each_user_channel_coefficient_H' * bs_filter_U(:,tmp_i));
                    elseif tmp_i == tmp_j
                        psi_matrix(tmp_i,tmp_j)  = 0;
                    end
                end
            end
    
            % equation (9) - geometric programming part
            gpvar variable_vector_p(total_data_num_N_d);
            %variable_vector_p = ones(total_data_num_N_d,1);
            tmp_value = (diagnal_matrix_D + scailed_factor_beta.^2 * psi_matrix) * variable_vector_p + sigma_n^2 * scailed_factor_beta.^2 * ones(total_data_num_N_d,1);
            problem = 1;
            for tmp_i = 1:total_data_num_N_d
                problem = problem * variable_vector_p(tmp_i)^(-1) * tmp_value(tmp_i);
            end
            constr = [variable_vector_p.^0.5' * variable_vector_p.^0.5<=P_max];
    
            [min_value,solution,status] = gpsolve(problem,constr,'min');
            assign(solution);
            power_allocation_P = diag(variable_vector_p);
    
            % equation (4), UL-receiver filter(DL-transmitter filter)
            bs_filter_U = [];
            scailed_factor_beta = [];
            for tmp_k = 1: User_Num_K
                U_k_beta_k = inv(each_user_channel_coefficient_H * user_filter_T * power_allocation_P * user_filter_T' * each_user_channel_coefficient_H' + sigma_n^2*eye(Tx_antenna_num_N_T)) * each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k)  *  user_filter_T(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k) * power_allocation_P(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k * Rx_antenna_num_N_R_k,1+(tmp_k-1) * (Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k); 
                beta_k = diag(sqrt(sum(abs(U_k_beta_k).^2)));
                U_k = U_k_beta_k*inv(beta_k);
                bs_filter_U = [bs_filter_U,U_k];
                scailed_factor_beta = blkdiag(scailed_factor_beta,beta_k);
            end 
        
            %% downlink channel
            % equation (24)
            mse_epsilon = diag(diag(diagnal_matrix_D) + diag(scailed_factor_beta.^2) ./ diag(power_allocation_P) .* (psi_matrix * diag(power_allocation_P)) + sigma_n^2 * (diag(scailed_factor_beta).^2) ./ (diag(power_allocation_P)));
            %tmp_i = 1;
            %tmp_2 = scailed_factor_beta(tmp_i,tmp_i)/power_allocation_P(tmp_i,tmp_i)* psi_matrix * diag(power_allocation_P);
            %diagnal_matrix_D(tmp_i, tmp_i) + tmp_2(tmp_i) + sigma_n^2 * scailed_factor_beta(tmp_i,tmp_i)^2 / (power_allocation_P(tmp_i,tmp_i))
            %tmp_value = (diagnal_matrix_D + scailed_factor_beta.^2 * psi_matrix) * diag(power_allocation_P) + sigma_n^2 * scailed_factor_beta.^2 * ones(total_data_num_N_d,1);
            %tmp_value./diag(power_allocation_P)
            
            
            % equation (11)
            power_allocation_Q = diag(sigma_n^2 * (mse_epsilon-diagnal_matrix_D - (scailed_factor_beta.^2)*transpose(psi_matrix))^-1 * scailed_factor_beta.^2 * ones(total_data_num_N_d,1) );
            
            % equation (10), DL-receiver filter (UL-transmitter filter)
            user_filter_T = [];
            scailed_factor_beta_tmp = [];
            for tmp_k = 1: User_Num_K
                T_k_beta_k_tmp = inv(each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k)' * bs_filter_U * power_allocation_Q * bs_filter_U' * each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k) + sigma_n^2*eye(Rx_antenna_num_N_R_k)) * each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k)'  *  bs_filter_U(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k) * power_allocation_Q(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k * Rx_antenna_num_N_R_k,1+(tmp_k-1) * (Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k); 
                beta_k_tmp = diag(sqrt(sum(abs(T_k_beta_k_tmp).^2)));
                T_k = T_k_beta_k_tmp*inv(beta_k_tmp);
                user_filter_T = blkdiag(user_filter_T,T_k); %[bs_filter_U,U_k];
                scailed_factor_beta_tmp = blkdiag(scailed_factor_beta_tmp,beta_k_tmp);
            end 
    
            %% uplink channel
            % equation (12)
            power_allocation_P = diag(sigma_n^2 * (mse_epsilon-diagnal_matrix_D - (scailed_factor_beta_tmp.^2)*transpose(psi_matrix))^-1 * scailed_factor_beta_tmp.^2 * ones(total_data_num_N_d,1) );
            
            % equation (4), UL-receiver filter(DL-transmitter filter)
            bs_filter_U = [];
            scailed_factor_beta = [];
            for tmp_k = 1: User_Num_K
                U_k_beta_k = inv(each_user_channel_coefficient_H * user_filter_T * power_allocation_P * user_filter_T' * each_user_channel_coefficient_H' + sigma_n^2*eye(Tx_antenna_num_N_T)) * each_user_channel_coefficient_H(:,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k)  *  user_filter_T(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k,1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k) * power_allocation_P(1+(tmp_k-1)*(Rx_antenna_num_N_R_k):tmp_k * Rx_antenna_num_N_R_k,1+(tmp_k-1) * (Rx_antenna_num_N_R_k):tmp_k*Rx_antenna_num_N_R_k); 
                beta_k = diag(sqrt(sum(abs(U_k_beta_k).^2)));
                U_k = U_k_beta_k*inv(beta_k);
                bs_filter_U = [bs_filter_U,U_k];
                scailed_factor_beta = blkdiag(scailed_factor_beta,beta_k);
            end 
            %diag(diagnal_matrix_D) + diag(scailed_factor_beta.^2) ./ diag(power_allocation_P) .* (psi_matrix * diag(power_allocation_P)) + sigma_n^2 * (diag(scailed_factor_beta).^2) ./ (diag(power_allocation_P))
            sum_rate = -log2(prod(diag(diagnal_matrix_D) + diag(scailed_factor_beta.^2) ./ diag(power_allocation_P) .* (psi_matrix * diag(power_allocation_P)) + sigma_n^2 * (diag(scailed_factor_beta).^2) ./ (diag(power_allocation_P))))
        end
    end
end