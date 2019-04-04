function MainAlgorithm(sensors)
load('sub2s.mat'); 
load('B.mat');
%sensors = [1 0 1 1 1 0 0 1 1 0 0 0 0 1 1];

f = find(diff([0,sensors,0]==1)); %index of every change
x = f(1:2:end-1);  % Start indices of 1s
y = f(2:2:end)-x;  % Consecutive 1s counts
[M_maxO,d] = max(y); %Index in x of the location of the max consecutive line of 1s

if M_maxO <2
    error('M too small. This array cannot function as an NSA')
end

first = x(1); 
last = find(sensors,1,'last'); %Because sensors is 1s and 0s, this will find the last 1
L_orig = last-first+1;
sensors = sensors(first:last); %trims any outside zeros

badcase = 0;
match = 0; %Determines whether the configuration matches a known one
while true
    for L = L_orig:-1:7
        disp('Im checking a new aperture')
        lookhere = find(cell2mat(sub2s(:,2))==L);
        while isempty(lookhere)==0 %there is a valid config of that L 
            for i = lookhere(1):lookhere(end)
                sub2s_check = sub2s{i,1};          
                for M_max = M_maxO:-1:2 % 2 is the lowest M we will allow
                    disp('Im checking a new M')
                    while M_max==(B(i,1)*(1+B(i,2)))%M*(1+p)
                    disp('Im testing different p values')    
                    PF = zeros(1,length(sub2s_check)); %pass/fail matrix
                        for j = 1:length(sub2s_check)
                           if sensors(j)==0 && sub2s_check(j)==1
                           PF(j)=1; %A sensor failure in a critical location
                           % sum(PF)=0 if it passes sensor by sensor
                           end       
                        end
                        if sum(PF)==0
                            match = 1;
                            M = B(i,1);
                            p = B(i,2);
                            l = B(i,3); %aperture
                            
                            break % out of M_max==B(i,1)*(1+B(i,2))%M*(1+p) 
                        end                   
                    end
                if match == 1
                  % out of for M_max = M_maxO:-1:2 
                    break
                end
                end
            if match == 1
              % out of i = lookhere(1):lookhere(end)
                break
            end
            end
        if match == 1
          % out of isempty(lookhere)==0
            break
        end
        end        
    if match == 1
    % out of L = L_orig:-1:7
        break
    end
    end
   
if match == 1
    % out of while true -> end function
    break
end             
%%%%%%%%%%%%%%%%%%%%%%%%% NOT AS GOOD AS A ULA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for L = L_orig:-1:7
        for M = 2:M_maxO %M
        sub2 = zeros(1,L);
        sub2(1:M:end) = 1; % Required sub 2 with a certain M and L
            if mod(M_maxO,M)==0
                p = M_maxO/M-1;
                PF = zeros(1,length(sub2)); %pass/fail matrix
                for j = 1:length(sub2s)
                   if sensors(j)==0 && sub2(j)==1
                   PF(j)=1; %A sensor failure in a critical location
                   % sum(PF)=0 if it passes sensor by sensor
                   end       
                end
            if sum(PF)==0
                break %out to L=L_orig:-1:7
            end
            end
        end    
    end
if badcase==1
    break % out of while true loop
end 
break %out of while true loop to disp('Array is not able to be analyzed as an NSA')
end %end of while true loop

%%%% Plot!!
if match == 1
    NSA_Analysis(1,0,2,1,M,p,l);
elseif badcase == 1
    NSA_Analysis(1,0,2,1,M,p,L);
    msgbox('Warning: This configuration is NOT guaranteed','','warn')
else
    disp('Array is not able to be analyzed as an NSA')
end
end
