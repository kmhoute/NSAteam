%Merged
% This script combines the changing SLH and testing all M and p. This
% function loops through M, calculates and loops through p, and then tests
% all values of L. 
Mtotest = 32; %maximum nummber of M
results = cell(Mtotest-1,1); 
A = cell(Mtotest-1,1);
for M = 2:Mtotest
    max_p = ((floor(64/M))-1);
    results{M-1} = cell(max_p,2);
    for k = 1:max_p+1
        p_actual = k-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        peak_compare = zeros(64, 3);
        while (M+M*p_actual+1)>=64
            %disp('Error1')
            %disp(num2str(p_actual))
            %disp(num2str(M))
            max_p=max_p-1;
            break
        end
    for L = (M+M*p_actual+1):64  % comment this line to run beampattern only; 
                      % uncomment for prod vs min analysis on # of sensors
        [Bmin, Bprod, N] = ProductMinBeampatternC(M,p_actual,L,0);

        if (p_actual+1)*M>=N*M
            disp('Error2')
        end
    
        Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end));
        locmin_1_min = find(Bmin_MinPos~=0,1);
        locmin_1_min = round(length(Bmin)/2)+locmin_1_min-1;
        Bprod_MinPos = islocalmin(Bprod(floor((length(Bprod)/2)):end));
        locmin_1_prod = find(Bprod_MinPos~=0,1);
        locmin_1_prod = round(length(Bmin)/2)+locmin_1_prod-1;
        peaks_Bmin  = max(Bmin(locmin_1_min:end));
        peaks_Bprod = max(Bprod(locmin_1_prod:end));

        peak_compare(L, 2) = peaks_Bmin;
        peak_compare(L, 3) = peaks_Bprod;
        peak_compare(L,1) = L;
    end  %%%%% comment this line and the following lines to run beampattern only; 
       %%%% uncomment for prod vs min analysis on # of sensors
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    results{M-1}{k,1} = peak_compare;
    no_ofindexes = find(results{M-1}{k,1}(:,2)<=-13);
    results{M-1}{k,2}((1:length(no_ofindexes)),1) = M; %M value
    results{M-1}{k,2}((1:length(no_ofindexes)),2) = p_actual; %p value
    results{M-1}{k,2}((1:length(no_ofindexes)),3) = results{M-1}{k,1}(no_ofindexes,1);%L
    end
    A{M-1}=cat(1,results{M-1}{:,2});
end
B = cat(1,A{:}); %These give all possibilities for a successful array