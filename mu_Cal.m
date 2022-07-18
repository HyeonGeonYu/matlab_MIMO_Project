function mu = mu_Cal(Lambda_matrix,W_matrix,p_0)
    
    
  for iter_ = 1: 5 % rank(H)
    rho_vector = diag( Lambda_matrix ) .* diag( W_matrix );

    mu = ( (sum( diag( Lambda_matrix^(-0.5) )  .* diag( W_matrix^(0.5) ) ))...
     / ( p_0 + sum( diag( Lambda_matrix^(-1) ) ) ) )^2;
    
    
    if all(rho_vector - mu>0)
      break
    end
    
  end

  if ~all(rho_vector - mu>0) 
     error("mu value does not exist")
  end
end