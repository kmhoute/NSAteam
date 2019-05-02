function MainAlgorithm_realdata(sensors) %, configurations)
%     if strcmp(configurations, 'old') == 1
%        % do nothing so no 
%     elseif strcmp(configurations, 'new') == 1
%         AllArrays(stuff) % put actual inputs
%     else
%         disp('Error: configurations must be "new" or "old"')
%     end
    
    load('sub2s.mat'); 
    load('finalConfigs_good.mat');
    load('sub2s_bad.mat'); 
    load('finalConfigs_bad.mat');
    
    %sensors = [1 0 1 1 1 0 0 1 1 0 0 0 0 1 1];
    %          [1 0 0 0 1 1 1 1 1 0 0 0 1 0 0 1 1 1 0] is a 'bad case'
    %          [1 1 1 1 1 1 0 0 0 0 1 0 0 1] is a case to check that it decreases L
    %          [1 0 1 1 1 0 0 1 1 0 0 0 0 1 1] is not an NSA


    f = find(diff([0,sensors,0]==1)); % index of every change
    x = f(1:2:end-1);  % Start indices of 1s
    y = f(2:2:end)-x;  % Consecutive 1s counts
    [maxM,~] = max(y); %Index in x of the location of the max consecutive line of 1s

    if maxM <2
        error('M is too small. This array cannot function as an NSA.')
    end

    first = x(1); 
    last = find(sensors,1,'last'); %Because sensors is 1s and 0s, this will find the last 1
    L_orig = last-first+1;
    sensors = sensors(first:last); %trims any outside zeros

    badcase = 0;
    match = 0; %Determines whether the configuration matches a known one
    while true
        for L = L_orig:-1:7  % check aperture decreasing
            disp('Im checking a new aperture, good'); disp(num2str(L))
            sub2aperL = find(cell2mat(sub2s(:,2))==L);  % saves indices of sub2s with aperture L
            if isempty(sub2aperL) == 0 % there is a valid config of that L
                lookhere = find(finalConfigs_good(sub2aperL(1):sub2aperL(end),4)<=maxM) + sub2aperL(1)-1; % with valid sub1+ext
                if isempty(lookhere)==0  % there is a valid config of that L (lookhere is not empty) 
                    for i = lookhere(1):lookhere(end)  % search through the valid sub2s 
                        sub2s_check = sub2s{i,1};  % saves current sub2 being checked 
                        PF = zeros(1,length(sub2s_check)); %pass/fail matrix
                        for j = 1:length(sub2s_check)
                            if sensors(j)==0 && sub2s_check(j)==1
                                PF(j)=1; %A sensor failure in a critical location
                           % sum(PF)=0 if it passes sensor by sensor
                            end       
                        end
                        if sum(PF)==0
                            match = 1;
                            M = finalConfigs_good(i,1);
                            p = finalConfigs_good(i,2);
                            disp(finalConfigs_good(i,:));
                            %l = finalConfigs_good(i,3); %aperture

                            break % out of i = lookhere(1):lookhere(end) 
                        end                   
                    end
                end
            end
            if match == 1
                break % out of L = L_orig:-1:7
            end
        end     
        if match == 1 
            break % out of while true -> end function
        end             
    %%%%%%%%%%%%%%%%%%%%%%%%% NOT AS GOOD AS A ULA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for L = L_orig:-1:7  % check aperture decreasing
            disp('Im checking a new aperture, bad'); disp(num2str(L))
            sub2aperL = find(cell2mat(sub2s_bad(:,2))==L);  % saves indices of sub2s with aperture L
            if isempty(sub2aperL) == 0 % there is a valid config of that L
                lookhere = find(finalConfigs_bad(sub2aperL(1):sub2aperL(end),4)<=maxM) + sub2aperL(1)-1; % with valid sub1+ext
                if isempty(lookhere)==0  % there is a valid config of that L (lookhere is not empty) 
                    for i = lookhere(1):lookhere(end)  % search through the valid sub2s 
                        sub2s_check = sub2s_bad{i,1};  % saves current sub2 being checked 
                        PF = zeros(1,length(sub2s_check)); %pass/fail matrix
                        for j = 1:length(sub2s_check)
                            if sensors(j)==0 && sub2s_check(j)==1
                                PF(j)=1; %A sensor failure in a critical location
                           % sum(PF)=0 if it passes sensor by sensor
                            end       
                        end
                        if sum(PF)==0
                            badcase = 1;
                            M = finalConfigs_bad(i,1);
                            p = finalConfigs_bad(i,2);
                            disp(finalConfigs_bad(i,:));
                            %l = finalConfigs_bad(i,3); %aperture

                            break % out of i = lookhere(1):lookhere(end) 
                        end                   
                    end
                end
            end
            if badcase == 1
                break % out of L = L_orig:-1:7
            end
        end      
        break % out of while true -> end function
    end %end of while true loop

    %%%% Plot!!
    if match == 1
        NSA_Analysis(0,1,2,1,M,p,L); %figure;
        BP_Formation.Nested(M,p,L,1);
    elseif badcase == 1
        NSA_Analysis(0,1,2,1,M,p,L); %figure;
        BP_Formation.Nested(M,p,L,1);
        msgbox('This configuration is NOT guaranteed','Warning','warn')
    else
        msgbox('Array is not able to be analyzed as an NSA', 'Error', 'error')
    end
end
