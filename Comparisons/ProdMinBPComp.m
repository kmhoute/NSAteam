function [u, B1, B2, Bmin, Bprod, N, add1] = ProdMinBPComp(M,p,L)
%%% M and N are the number of sensors in the basic subarrays
%%% add1 and add2 are the additional sensors in subarray1 and subarray2
%%% U1 and U2 are undersampling factors

    N = ceil(L/M); % calculates the number of sensors in N
    U1 = 1; % this will not change
    U2 = M; % will always be M since we now use add1 to adjust from basic nested
    add1 = M*p; % adds sensors to subarray1 in multiples of M
    add2 = 0; % will always be zero

%%% the rest uses the ProductMinBeamPattern function code without the plotting %%%    
    uDelta = 0.001;
    u = -1:uDelta:1;
    
    B1 = BP_Formation.Linear(0, M+add1, 0, U1*0.5);
    B1 = B1/max(abs(B1));
    
    
    B2 = BP_Formation.Linear(0, N+add2, 0, U2*0.5);
    B2 = B2/max(abs(B2));
   

    Bmin = min(abs([B1;B2]));
    Bmin = 20*log10(abs(Bmin));

    Bprod = B1.*conj(B2);
    Bprod = 10*log10(abs(Bprod));
    
    B1 = 20*log10(abs(B1));
    B2 = 20*log10(abs(B2));
 end
