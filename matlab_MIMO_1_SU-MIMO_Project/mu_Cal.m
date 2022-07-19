function [mu,diag_f_matrix] = mu_Cal(Lambda_matrix,W_matrix,p_0,B)
    
    
  for iter_ = 0: B-1 %
    
    Lambda_vector = diag( Lambda_matrix);
    W_vector =diag( W_matrix );
    rho_vector = Lambda_vector .* W_vector;

    mu = ( (sum( Lambda_vector(1:B-iter_).^(-0.5) .* W_vector(1:B-iter_).^(0.5)  ))...
     / ( p_0 + sum( Lambda_vector(1:B-iter_).^(-1) ) ) )^2;
    
    
    if all(rho_vector(1:B-iter_)  - mu>0)
      break
    end
    
  end

  diag_f_matrix = ( mu^(-0.5) * Lambda_matrix^(-0.5) * W_matrix^(0.5)...
   - (Lambda_matrix)^(-1))^(0.5);
  
  diag_f_tmp = diag(diag_f_matrix);
  diag_f_tmp(B+1-iter_:B) = 0 ;
  diag_f_matrix = diag(diag_f_tmp);

  diag_f_matrix(diag_f_matrix<0) = 0;

  if ~all(rho_vector(1:B-iter_)  - mu>0)
     error("mu value does not exist")
  end
end