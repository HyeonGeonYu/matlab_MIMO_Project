function [H] = Channel_Gen_ML(SP)

    Nr = SP.Nr;
    Nu = SP.Nu;
    H_type = SP.H_type;
    L = SP.L;
    switch H_type
        
        case 'Rayleigh'
            H = 1/sqrt(2)*(randn(Nr,Nu) + 1j*randn(Nr,Nu));
    end
    
    end