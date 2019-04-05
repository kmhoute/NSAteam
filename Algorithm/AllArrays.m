% AllArrays
% This script combines the changing SLH and testingallMandp. This
% function loops through M, calculates and loops through p, and then tests
% all values of L. 
Mtotest = 63; %maximum nummber of M
results = cell(Mtotest-1,1); 
A = cell(Mtotest-1,1);
A_bad = cell(Mtotest-1,1);
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
        [~,~,~,Bmin, Bprod, N,~] = BP_Formation.Nested(M,p_actual,L,0);

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
        peak_compare(L,1) = M*(N-1)+1;
    end  %%%%% comment this line and the following lines to run beampattern only; 
       %%%% uncomment for prod vs min analysis on # of sensors
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    results{M-1}{k,1} = peak_compare;
    no_ofindexes = find(results{M-1}{k,1}(:,2)<=-13);
    results{M-1}{k,2}((1:length(no_ofindexes)),1) = M; %M value
    results{M-1}{k,2}((1:length(no_ofindexes)),2) = p_actual; %p value
    results{M-1}{k,2}((1:length(no_ofindexes)),3) = results{M-1}{k,1}(no_ofindexes,1);
    results{M-1}{k,2}((1:length(no_ofindexes)),4) = M*(p_actual+1);
    
    no_ofindexes_bad = find(results{M-1}{k,1}(:,2)>-13);
    results{M-1}{k,3}((1:length(no_ofindexes_bad)),1) = M; %M value
    results{M-1}{k,3}((1:length(no_ofindexes_bad)),2) = p_actual; %p value
    results{M-1}{k,3}((1:length(no_ofindexes_bad)),3) = results{M-1}{k,1}(no_ofindexes_bad,1);
    results{M-1}{k,3}((1:length(no_ofindexes_bad)),4) = M*(p_actual+1);
    end
    A{M-1}=cat(1,results{M-1}{:,2}); % save the good
    A_bad{M-1}=cat(1,results{M-1}{:,3}); % save the bad
end
%%%% good cases file %%%%
B = cat(1,A{:}); %These give all possibilities for a successful array
B = unique(B,'rows');
B = sortrows(B,[-3,-1,-2]); 
save('B.mat','B') %Saves the matrix 

sub2s = cell(length(B),2);
for i=1:length(B)
    sub2s{i,1} = false(1,B(i,3));
    sub2s{i,1}(1:B(i,1):end)=true;
    sub2s{i,2} = B(i,3);
end
save('sub2s.mat','sub2s') %Saves the matrix 

%%%% bad cases file %%%%
B_bad = cat(1,A_bad{:}); %These give all possibilities for a successful array
B_bad = unique(B_bad,'rows');
B_bad = sortrows(B_bad,[-3,-1,-2]); 
save('B_bad.mat','B_bad') %Saves the matrix 

sub2s_bad = cell(length(B_bad),2);
for i=1:length(B_bad)
    sub2s_bad{i,1} = false(1,B_bad(i,3));
    sub2s_bad{i,1}(1:B_bad(i,1):end)=true;
    sub2s_bad{i,2} = B_bad(i,3);
end
save('sub2s_bad.mat','sub2s_bad') %Saves the matrix 
