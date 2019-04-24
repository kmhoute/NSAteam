function [finalConfigs_good, sub2s, finalConfigs_bad, sub2s_bad] = AllArrays(maxL)
%%% This function combines the changing SLH and testingallMandp. 
%%%      maxL is the maximum aperture you want to generate
%%% This function loops through all M, loops through all L, and then calculates 
%%% and loops through p. 
%clear all; clc;
tic;
%maxM = 63;
%maxL = 64; 
    finalConfigs_good = zeros(1,4); % initiates the array (will become a matrix)
    finalConfigs_bad = zeros(1,4);  % initiates the array (will become a matrix)
    
    for M = 2:(maxL-1) % shortest M is 2
        for L = (M+1):maxL % shortest L is 1 more than M
            max_p = ((floor(L/M))-1); 
            if (M+M*max_p)>(L-1) % failsafe: makes sure sub1 is not equal to L
                %disp('Error1')
                %disp(num2str(p))
                %disp(num2str(M))
                %disp(num2str(L))
                max_p=max_p-1; % decrease p
            end
            
            for p = 0:max_p % loop through extension factor (indices are p+1 to start at 1 not 0)
                peak_compare = zeros((max_p+1), 2); % matrix for psl height check
                [u,~,~,Bmin, ~, N,~] = BP_Formation.Nested(M,p,L,0); % gennerates BP for the case
                
                % locates first null of the minimum processed B
                Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end)); % logic matrix
                locmin_1_min = find(Bmin_MinPos~=0,1); % returns the null's index in Bmin_MinPos
                locmin_1_min = round(length(Bmin)/2)+locmin_1_min-1; % returns the null's index in Bmin 
                
                peak_compare((p+1),1) = p; % saves the extension factor
                peak_compare((p+1),2) = max(Bmin(locmin_1_min:end)); % saves the max value which is the PSL height
            end % ends ext factor (p) loop  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            % this chunk stores configs for the different p values (with the same M and L)
            no_ofindexes = find(peak_compare(:,2)<=-13); % saves indices of successful configs
            Configs_good(1:length(no_ofindexes),1) = M; 
            Configs_good(1:length(no_ofindexes),2) = peak_compare(no_ofindexes,1);
            Configs_good(1:length(no_ofindexes),3) = L;
            Configs_good(1:length(no_ofindexes),4) = M*((peak_compare(no_ofindexes,1))+1); % stores totalM length
            
            % concatenates the cases stored above to the previously tested success cases
            finalConfigs_good = cat(1,finalConfigs_good(:,:), Configs_good(:,:));
           
            no_ofindexes_bad = find(peak_compare(:,2)>-13);
            Configs_bad(1:length(no_ofindexes_bad),1) = M;
            Configs_bad(1:length(no_ofindexes_bad),2) = peak_compare(no_ofindexes_bad,1);
            Configs_bad(1:length(no_ofindexes_bad),3) = L;
            Configs_bad(1:length(no_ofindexes_bad),4) = M*((peak_compare(no_ofindexes_bad,1))+1);
            
            finalConfigs_bad = cat(1,finalConfigs_bad(:,:),Configs_bad(:,:));
        end % ends L loop
    end % ends M loop
    
   %%%%% good cases file generation %%%%%
    finalConfigs_good = unique(finalConfigs_good,'rows'); % deletes any duplicate rows
    finalConfigs_good = sortrows(finalConfigs_good,[-3,-1,-2]); % sorts by L, M, then p
    finalConfigs_good = finalConfigs_good(1:(end-1),:); % deletes the zeros row from initation
    save('finalConfigs_good.mat','finalConfigs_good') % Saves the matrix to the file

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
