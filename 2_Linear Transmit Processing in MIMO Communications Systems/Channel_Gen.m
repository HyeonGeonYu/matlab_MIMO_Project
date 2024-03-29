function [H] = Channel_Gen(SP)

    Nr = SP.Nr;
    Nt = SP.Nt;
    H_type = SP.H_type;

    switch H_type
        
        case 'Rayleigh'
            H = 1/(2^1.5)*(randn(Nr,Nt) + 1j*randn(Nr,Nt));
    end
    
    end