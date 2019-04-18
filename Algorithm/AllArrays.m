function [finalConfigs_good, sub2s, finalConfigs_bad, sub2s_bad] = AllArrays(maxL)
%%% This function combines the changing SLH and testingallMandp. 
%%%      maxL is the maximum aperture you want to generate
%%% This function loops through all M, loops through all L, and then calculates 
%%% and loops through p. 
%clear all; clc;
tic;
%maxM = 63;
%maxL = 64; 
    finalConfigs_good = zeros(1,4);
    finalConfigs_bad = zeros(1,4);
    
    for M = 2:(maxL-1)
        for L = (M+1):maxL
            max_p = ((floor(L/M))-1);
            if (M+M*max_p)>(L-1)
                %disp('Error1')
                %disp(num2str(p))
                %disp(num2str(M))
                %disp(num2str(L))
                max_p=max_p-1;
            end
            
            for p = 0:max_p
                peak_compare = zeros((max_p+1), 2);
                [~,~,~,Bmin, ~, N,~] = BP_Formation.Nested(M,p,L,0);

                Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end));
                locmin_1_min = find(Bmin_MinPos~=0,1);
                locmin_1_min = round(length(Bmin)/2)+locmin_1_min-1;
                
                peak_compare((p+1),2) = max(Bmin(locmin_1_min:end));
                peak_compare((p+1),1) = p;
            end % ends ext factor (p) loop  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            no_ofindexes = find(peak_compare(:,2)<=-13);
            Configs_good(1:length(no_ofindexes),1) = M;
            Configs_good(1:length(no_ofindexes),2) = p;
            Configs_good(1:length(no_ofindexes),3) = L;
            Configs_good(1:length(no_ofindexes),4) = M*(p+1);
            
            finalConfigs_good = cat(1,finalConfigs_good(:,:), Configs_good(:,:));
           
            no_ofindexes_bad = find(peak_compare(:,2)>-13);
            Configs_bad(1:length(no_ofindexes_bad),1) = M;
            Configs_bad(1:length(no_ofindexes_bad),2) = p;
            Configs_bad(1:length(no_ofindexes_bad),3) = L;
            Configs_bad(1:length(no_ofindexes_bad),4) = M*(p+1);
            
            finalConfigs_bad = cat(1,finalConfigs_bad(:,:),Configs_bad(:,:));
        end % ends L loop
        
    end % ends M loop
   %%%%% good cases file %%%%%
    finalConfigs_good = unique(finalConfigs_good,'rows');
    finalConfigs_good = sortrows(finalConfigs_good,[-3,-1,-2]);
    finalConfigs_good = finalConfigs_good(1:(end-1),:);
    save('finalConfigs_good.mat','finalConfigs_good') %Saves the matrix 

    sub2s = cell(length(finalConfigs_good),2);
    for i=1:length(finalConfigs_good)
        sub2s{i,1} = false(1,finalConfigs_good(i,3));
        sub2s{i,1}(1:finalConfigs_good(i,1):end)=true;
        sub2s{i,2} = finalConfigs_good(i,3);
    end
    save('sub2s.mat','sub2s') %Saves the matrix 

   %%%%% bad cases file %%%%%
    finalConfigs_bad = unique(finalConfigs_bad,'rows');
    finalConfigs_bad = sortrows(finalConfigs_bad,[-3,-1,-2]);
    finalConfigs_bad = finalConfigs_bad(1:(end-1),:);
    save('finalConfigs_bad.mat','finalConfigs_bad') %Saves the matrix 

    sub2s_bad = cell(length(finalConfigs_bad),2);
    for i=1:length(finalConfigs_bad)
        sub2s_bad{i,1} = false(1,finalConfigs_bad(i,3));
        sub2s_bad{i,1}(1:finalConfigs_bad(i,1):end)=true;
        sub2s_bad{i,2} = finalConfigs_bad(i,3);
    end
    save('sub2s_bad.mat','sub2s_bad') %Saves the matrix 
    toc;
end
