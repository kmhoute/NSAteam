function MainAlgorithm(sensors)
load('sub2s.mat'); 
load('B.mat');
%sensors = [1 0 1 1 1 0 0 1 1 0 0 0 0 1 1];

f = find(diff([0,sensors,0]==1)); %index of every change
x = f(1:2:end-1);  % Start indices of 1s
y = f(2:2:end)-x;  % Consecutive 1s counts
[M_max,~] = max(y); %Index in x of the location of the max consecutive line of 1s
% position = x(i); If i need this, replace ~ with i
if M_max <2
    error('M too small. This array cannot function as an NSA')
end
Mdecmax = M_max-2;

first = x(1); 
last = find(sensors,1,'last');
Lorig = last-first+1;
sensors = sensors(first:last); %trims any outside zeros

plots = 0; %will use later to turn on the plots
working = 0; %will determine whether we are using a verified config (in B matrix)
match = 0; %
for k=0:(Lorig-7)
    L=Lorig-k;
    lookhere = find(cell2mat(sub2s(:,2))==L);
    while isempty(lookhere)==0 %there is a valid config of that L 
        for i = lookhere(1):lookhere(end)
            sub2s_check = sub2s{i,1};
            
            for m=1:(Mdecmax) % 2 is the lowest M we will allow
                if (B(i,1)==M_max) && (B(i,2)==0) 
                    M = M_max;
                    p = 0;
                    plots = 1; %Shows the plot
                    match = 1;
                elseif (B(i,2)~=0) && (mod(M_max,(1+B(i,2)))==0) 
                    M = floor(M_max/(1+B(i,2)));
                    p = B(i,2);
                    plots = 1; 
                    match = 1;
                else 
                    M_max = M_max-m;
                end
                
                if match == 1
                    break
                end
            end
                       
            PF = zeros(1,length(sub2s_check)); %pass/fail matrix
           for j = 1:length(sub2s_check)
               if sensors(j)==0 && sub2s_check(j)==1
               PF(j)=1; %A sensor failure in a critical location
               end       
           end

           if sum(PF)==0 && match==1
                NSA_Analysis(1,0,2,plots,M,p,L);
                %data: 1=simulated 0=real
                %filters: 1=applies a filter 0=does not apply a filter
                %processor: 1=product 2=minimum
                %plots: 1:plots the DOA and freq 0:does not plot
                %L
                disp(M); disp(p); disp(L); disp(i);
                
            elseif (plots == 0) || (match == 0)
               %plot the result knowing that it will not be a passing
               %config. Assume no extension
               plots=1;
               NSA_Analysis(1,0,2,plots,M_max,0,L);
               msgbox('Warning: This configuration is NOT guaranteed')
            end     
%            elseif L==7 && sum(PF)>0
%                msgbox('No Configurations were found')
        break
        end
        break
    end 
    break
end
end