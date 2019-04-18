%function AllArrays(maxM,maxL)
%%% This function combines the changing SLH and testingallMandp. 
%%%   maxM is the maximum length you want to test
%%%        typically this will be maxL-1
%%%   maxL is the maximum aperture you want to generate
%%% This function loops through all M, calculates and loops through p, and then tests
%%% all values up to maxL. 
tic;
maxM = 63;
maxL = 64;
    % results = cell(maxM-1,1); 
    finalConfigs_good = zeros(1,4);
    finalConfigs_bad = zeros(1,4);
    for M = 2:maxM
        max_p = ((floor(maxL/M))-1);
        extFactor = [max_p,2]; %results{M-1} = cell(max_p,2);
        for k = 1:max_p+1
            p_actual = k-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            peak_compare = zeros(maxL, 2);
            while (M+M*p_actual+1)>maxL
                %disp('Error1')
                %disp(num2str(p_actual))
                %disp(num2str(M))
                max_p=max_p-1;
                break
            end
            for L = (M+M*p_actual+1):maxL  
                [~,~,~,Bmin, ~, N,~] = BP_Formation.Nested(M,p_actual,L,0);

                if (p_actual+1)*M>=N*M
                    disp('Error2')
                end

                Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end));
                locmin_1_min = find(Bmin_MinPos~=0,1);
                locmin_1_min = round(length(Bmin)/2)+locmin_1_min-1;
                
                peak_compare(L,2) = max(Bmin(locmin_1_min:end));
                peak_compare(L,1) = L;
            end  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            no_ofindexes = find(peak_compare(:,2)<=-13);
            Configs_good(1:length(no_ofindexes),1) = M;
            Configs_good(1:length(no_ofindexes),2) = p_actual;
            Configs_good(1:length(no_ofindexes),3) = peak_compare(no_ofindexes,1);
            Configs_good(1:length(no_ofindexes),4) = M*(p_actual+1);
            
            finalConfigs_good = cat(1,finalConfigs_good(:,:), Configs_good(:,:));
            %results{M-1}{k,1} = peak_compare;
            %no_ofindexes = find(results{M-1}{k,1}(:,2)<=-13);
            %results{M-1}{k,2}((1:length(no_ofindexes)),1) = M; %M value
            %results{M-1}{k,2}((1:length(no_ofindexes)),2) = p_actual; %p value
            % results{M-1}{k,2}((1:length(no_ofindexes)),3) = results{M-1}{k,1}(no_ofindexes,1);
            %results{M-1}{k,2}((1:length(no_ofindexes)),4) = M*(p_actual+1);

            no_ofindexes_bad = find(peak_compare(:,2)>-13);
            Configs_bad(1:length(no_ofindexes_bad),1) = M;
            Configs_bad(1:length(no_ofindexes_bad),2) = p_actual;
            Configs_bad(1:length(no_ofindexes_bad),3) = peak_compare(no_ofindexes_bad,1);
            Configs_bad(1:length(no_ofindexes_bad),4) = M*(p_actual+1);
            
            finalConfigs_bad = cat(1,finalConfigs_bad(:,:),Configs_bad(:,:));
            %no_ofindexes_bad = find(results{M-1}{k,1}(:,2)>-13);
            %results{M-1}{k,3}((1:length(no_ofindexes_bad)),1) = M; %M value
            %results{M-1}{k,3}((1:length(no_ofindexes_bad)),2) = p_actual; %p value
            %results{M-1}{k,3}((1:length(no_ofindexes_bad)),3) = results{M-1}{k,1}(no_ofindexes_bad,1);
            %results{M-1}{k,3}((1:length(no_ofindexes_bad)),4) = M*(p_actual+1);
        end
        %tempConfigs_good{M-1}=cat(1,results{M-1}{:,2}); % save the good
        %tempConfigs_bad{M-1}=cat(1,results{M-1}{:,3}); % save the bad
    end
   %%%%% good cases file %%%%%
    %Configs_good = cat(1,tempConfigs_good{:}); %These give all possibilities for a successful array
    finalConfigs_good = unique(finalConfigs_good,'rows');
    finalConfigs_good = sortrows(finalConfigs_good,[-3,-1,-2]); 
    save('finalConfigs_good.mat','finalConfigs_good') %Saves the matrix 

    sub2s = cell(length(finalConfigs_good),2);
    for i=1:length(finalConfigs_good)
        sub2s{i,1} = false(1,finalConfigs_good(i,3));
        sub2s{i,1}(1:finalConfigs_good(i,1):end)=true;
        sub2s{i,2} = finalConfigs_good(i,3);
    end
    save('sub2s.mat','sub2s') %Saves the matrix 

   %%%%% bad cases file %%%%%
    %Configs_bad = cat(1,tempConfigs_bad{:}); %These give all possibilities for a successful array
    finalConfigs_bad = unique(finalConfigs_bad,'rows');
    finalConfigs_bad = sortrows(finalConfigs_bad,[-3,-1,-2]); 
    save('finalConfigs_bad.mat','finalConfigs_bad') %Saves the matrix 

    sub2s_bad = cell(length(finalConfigs_bad),2);
    for i=1:length(finalConfigs_bad)
        sub2s_bad{i,1} = false(1,finalConfigs_bad(i,3));
        sub2s_bad{i,1}(1:finalConfigs_bad(i,1):end)=true;
        sub2s_bad{i,2} = finalConfigs_bad(i,3);
    end
    save('sub2s_bad.mat','sub2s_bad') %Saves the matrix 
    toc;
%end