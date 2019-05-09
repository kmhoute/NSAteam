classdef ArraySystem
%%%% This class contains functions for forming beampatterns

%%%%%% The following functions produce plots for analysis
%%%%%%%%%% the 'TemporalSpatial' function creates the Temporal and Spatial Spectrum estimates
%%%%%%%%%% the 'NestedBP' function creates the Nested Sensor Array Beampattern
%%%%%%%%%% the 'UniformBP' function creates the Uniform Linear Array Beampattern
%%%%%%%%%% the 'LinearBP' function creates the Linear function to use for forming beampattern

%%%%%% The following functions generate configurations for the sparse array
%%%%%%%%%% the 'GenerateConfigurations' function generates the NSA configurations, good and bad
%%%%%%%%%% the 'SelectConfiguration' function analyzes configurations and selects best option

    methods(Static)
    %%%%%%% The following functions produce plots for analysis
        % creates the Temporal and Spatial Spectrum estimates
        function TemporalSpatial(data,filters,processor,plots,M,p,L,SNRdB,angle,f)
        %%%% This function uses:   
        %%%%%%%      real    data if data = 0  -OR-  simulated data if data = 1
        %%%% AND can apply a filter to the array:
        %%%%%%%      no filter  --> filters = 0  -OR-  with filter --> filters = 1
        %%%% AND applies a processor to the array:
        %%%%%%%  NO processor (ULA) --> processor = 0  -OR-
        %%%%%%%      product --> processor = 1  -OR-   minimum --> processor = 2
        %%%% AND applies an option for plots
        %%%%%%%      no plots --> plots = 0  -OR-  yes plots --> plots = 1 
        %%% ******* simulation data parameters can be ******* 
        %%%   *** added to input arguments if necessary ***
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%% data %%%%%%%%%%%%%%%
            if (data == 0)
                %%%% load path is dependant on the computer being used
                load %%%%%%%%'C:\Users\MaxSlezak\Google Drive\Classes\2018-2019\Senior Design\Measurement_Data\3\138'
                [totalData order]=Rearrange_Data(totalData);
                totalData=totalData(100:40100,3:L);
            elseif (data == 1)
                %%%%%%% Some constant parameters
                %angle = -34;
                %f = 8923.8;%%%%Signal frequency in Hz
                c = 343;%%%%%Signal speed in m/s
                deltat = 1/44100;%%%%Temporal sampling interval in s
                %SNRdB = 5;%%%%SNR in dB
                vars = 1;%%%%Signal variance
                varn = vars*10^(-SNRdB/10);%%%%Noise variance

                %%%%To create x = (exp(1j*2*pi*f*t-1j*pi*cosd(-34)*indices)+ ...
                %%%%exp(1j*2*pi*(-f)*t-1j*pi*cosd(-34)*indices)), make a matrix with the t
                %%%%values and indices values first
                times = (0:deltat:1)';
                locations = 0:(L-1);
                [indices,t] = meshgrid(locations,times);
                x = exp(1j*(2*pi*f*t-pi*cosd(angle)*indices))+conj(exp(1j*(2*pi*f*t-pi*cosd(angle)*indices)));
                clear indices t times;
                %%%%To the matrix x above, we need to add white noise
                totalData = x + sqrt(varn/2)*randn(size(x)) + 1i*sqrt(varn/2)*randn(size(x));
            else
                disp('Error: data = 0 or 1')
            end

            %%%%%%%%%%%%%%%%% filter %%%%%%%%%%%%%%%%%%%
            if (filters == 1)
                %Apply a highpass filter to remove noise
                band = designfilt('highpassiir','FilterOrder',50,...
                'PassbandFrequency', 3500, 'PassbandRipple',.2,...
                'SampleRate',44100); %This filter passes frequencies above 4000
                totalData=filter(band,totalData);
            elseif (filters == 0)
            else
                disp('Error: filter = 0 or 1')
            end

            %%%%%%%%%%%%%%% processor %%%%%%%%%%%%%%%
            if processor == 0
                F = fft(totalData,(L-2)*10,2);
                mF = mean(abs(F));
                mF = 20*log10(fliplr(fftshift(mF)/max(abs(mF))));
                w = linspace(-1,1,length(mF));
            elseif (processor == 1)|(processor == 2)
                sparseindices1 = (1:(M+M*p));%%% We add 1 at the end because MATLAB
                                            %%% index starts at 1, not 0
                sparseindices2 = (1:M:L);%%% We add 1 at the end because MATLAB
                                                 %%% index starts at 1, not 0
                data1 = zeros(size(totalData));
                data2 = zeros(size(totalData));
                data1(:,sparseindices1) = totalData(:,sparseindices1);
                data2(:,sparseindices2) = totalData(:,sparseindices2);
            
                if (processor==1)
                    %%%%Find the product
                    F = fft(data1,L*10,2) .* conj(fft(data2,L*10,2));
                    %%%%Find the average
                    mF = mean(abs(F));
                    %%%%Normalize, fftshift, flip left to right and convert to dB
                    mF = 10*log10(fliplr(fftshift(mF)/max(abs(mF))));
                elseif (processor==2)
                    %%%%Find the minimum of the absolute values
                    F = min(abs(fft(data1,L*10,2)) , abs(fft(data2,L*10,2)));
                    %%%%Find the average
                    mF = mean(abs(F));
                    %%%%Normalize, fftshift, flip left to right and convert to dB
                    mF = 20*log10(fliplr(fftshift(mF)/max(abs(mF))));
                end
                w = linspace(-1,1,length(mF));
            else 
                disp('Error: Input argument must be 1 or 2.');
                return;
            end

            %%%%%%%%%%%%%%% analysis plots %%%%%%%%%%%%%%%
            if (plots == 1)
                %%%%%%%%%%% Plot the spatial spectrum %%%%%%%%%%%%
                figure;
                subplot(2,1,1);
                plot(acosd(w),mF,'LineWidth',2);
                set(gca,'FontSize',12);
                ylabel('Power dB','FontWeight','bold', 'FontSize',16);
                xlabel('\theta centered','FontWeight','bold', 'FontSize',16);
                title('Spatial Spectral Estimation','FontWeight','bold','FontSize',20);
                xlim([0 180]); ylim([-20 0]); grid('on');

                %%%%%%%%%%% Plot the temporal spectrum %%%%%%%%%%%%
                F = fft(totalData,size(totalData,1)*10);
                mF = mean(abs(F),2);
                mF = fftshift(mF)/max(abs(mF));
                w = linspace(-1,1,length(mF));

                hold on; subplot(2,1,2);
                plot(w,20*log10(mF),'LineWidth',2);
                set(gca,'FontSize',12);
                ylabel('Power dB','FontWeight','bold', 'FontSize',16);
                xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold', 'FontSize',16);
                title('Temporal Spectral Estimation','FontWeight','bold', 'FontSize',20);
                xlim([-1 1]); ylim([-40 0]); grid('on');

            elseif (plots == 0)
            else
                disp('Error; plots must be 0 or 1')
            end
        end
    
        % creates the Nested Sensor Array Beampattern
        function [u, B1, B2, Bmin, Bprod, N, add1] = NestedBP(M, p, L, plots)
        %%% Nested creates a Nested Sensor Array Beampattern with 
        %%% M elements in subarray 1, N elements in subarray 2
        %%% p is the extension factor (applied to subarray 1)
        %%%       subarray total elements = M + M*p
        %%% L is the overarching aperture (range) of possible elements
        %%% set plots = 1 the beampattern plot, plots = 0 for no plot
            %%% M and N are the number of sensors in the basic subarrays
            M = M;
            N = ceil(L/M);  
        %%% add1 and add2 are the additional sensors in subarray1 and subarray2
            add1 = M*p;
            add2 = 0;
        %%% U1 and U2 are undersampling factors
            U1 = 1;
            U2 = M;

            uDelta = 0.001;
            u = -1:uDelta:1;
            %%%% LinearBP arguments: 
            %%%% 1. 0 for no plots %%% 2. elements in subarray 1 
            %%%% 3. extra sensors in subarray %%% 4. sensor spacing (lambda=0.5)
            [B1, ~] = ArraySystem.LinearBP(0, M+add1, 0, U1*0.5);
             B1 = B1/max(abs(B1));
            [B2, ~] = ArraySystem.LinearBP(0, N+add2, 0, U2*0.5);
             B2 = B2/max(abs(B2));

            Bmin = min(abs([B1;B2]));
            Bmin = 20*log10(abs(Bmin));

            Bprod = B1.*conj(B2);
            Bprod = 10*log10(abs(Bprod));
                     
            if plots == 1
                figure; 
                set(gca, 'FontSize', 12);
                plot(u, 20*log10(abs(B1)), 'LineStyle', '-.', 'LineWidth', 1, 'Color', 'Blue');
                hold on;
                plot(u, 20*log10(abs(B2)), 'LineStyle', '-.', 'LineWidth', 3, 'Color', 'Red');
                hold on; grid on;
                plot(u, Bmin, 'LineWidth', 2, 'Color', 'Black');
                hold on;
                plot(u, Bprod, '--','LineWidth', 1, 'Color', [0 0.6 0]);
                hold on;
                plot([-1 1],[-13 -13],'--','LineWidth',1,'Color','k');
                xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
                ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
                legend('Subarray 1','Subarray 2', 'Min', 'Product');
                xlim([-1 1]);
                ylim([-30 0]);
                title(sprintf('NSA Beampattern with L = %d, M = %d, add1 = %d ',L,M,add1), 'FontSize', 20, 'FontWeight', 'Bold')
            elseif plots == 0
            else
                error('plots must be 0 or 1')
            end
        end
        
        % creates the Uniform Linear Array Beampattern
        function [u, B, L] = UniformBP(L, plots)
        %%% UniformBP creates a Uniform Linear Array Beampattern with 
        %%% L is the overarching aperture (range) of possible elements
        %%% set plots = 1 the beampattern plot, plots = 0 for no plot

            uDelta = 0.001;
            u = -1:uDelta:1;
            %%%% LinearBP arguments: 
            %%%% 1. 0 for no plots %%% 2. elements in subarray 1 
            %%%% 3. extra sensors in subarray %%% 4. sensor spacing (lambda=0.5)
            [B, ~] = ArraySystem.LinearBP(0, L, 0, 0.5);
             B = B/max(abs(B));
             B = 20*log10(abs(B));
                     
            if plots == 1
                figure; grid on; 
                set(gca, 'FontSize', 12);
                plot(u, B, 'LineWidth', 2, 'Color', 'Black');
                plot([-1 1],[-13 -13],'--','LineWidth',1,'Color','k');
                xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
                ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
                xlim([-1 1]);
                ylim([-30 0]);
                title(sprintf('ULA Beampattern with Aperture = ',L))
            elseif plots == 0
            else
                error('plots must be 0 or 1')
            end
        end
        
        % creates the Linear function to use for forming beampattern
        function [B, u] = LinearBP(flag, N, varargin)
        %%%If flag is 1, plot the graphs. 
        %%%varargin: 1. shift = 0 for broadside steered, else enter u you want the
        %%%array to steer to 2. intersensor spacing in terms of wavelength 3.
        %%%weights: If you want to use matlab windows, use 'periodic' flag
        %%%whenever available Example: w = hann(N, 'periodic').' and use w while
        %%%calling the function 4. 5. xlimits
            m = length(varargin);
            OptArgs = {0 0.5  -1 1 1/N*ones(1,N)};
            OptArgs(1:m) = varargin;
            [shift, d, lu, uu, w] = OptArgs{:};
            n = (0:N-1)';
            u = lu:0.001:uu;
            B = zeros(size(u));
            for idx = 1:length(u)
                v = exp(1i*(n-(N-1)/2)*2*pi*(u(idx)-shift)*d);
                B(idx) = conj(w)*v;
            end
            % Normalize by the absolute value at u = 0
            B = B/sum(w);
            if flag
                %f = figure;
                set(gca, 'FontSize', 12);
                plot(u, (real(B)), 'LineWidth', 3, 'Color', 'Black');
                hold on;
                plot(u, (imag(B)), 'LineWidth', 3, 'Color', 'red');
                title('Beam Pattern in U - Space', 'FontWeight', 'Bold','FontSize', 20);
                grid on;
                xlabel('u', 'FontSize', 16, 'FontWeight', 'Bold');
                ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
                xlim([-1 1]);
            end
        end
        
    %%%%%%% The following functions generate configurations for the sparse array    
        % generates the NSA configurations, good and bad
        function [finalConfigs_good, sub2s, finalConfigs_bad, sub2s_bad] = GenerateConfigurations(maxL)
        %%% This function combines the changing SLH and testingallMandp. 
        %%%      maxL is the maximum aperture you want to generate
        %%% This function loops through all M, loops through all L, and then calculates 
        %%% and loops through p. 
            finalConfigs_good = zeros(1,4); % initiates the array (will become a matrix)
            finalConfigs_bad = zeros(1,4);  % initiates the array (will become a matrix)

            for M = 2:(maxL-1) % shortest M is 2
                for L = (M+1):maxL % shortest L is 1 more than M
                    max_p = ((floor(L/M))-1); 
                    if (M+M*max_p)>(L-1) % failsafe: makes sure sub1 is not equal to L
                        max_p=max_p-1; % decrease p
                    end
                    peak_compare = zeros((max_p+1), 3); % matrix for psl height check
                    for p = 0:max_p % loop through extension factor (indices are p+1 to start at 1 not 0)
                        [~,~,~,Bmin, ~, N,~] = ArraySystem.NestedBP(M,p,L,0); % gennerates BP for the case

                        % locates first null of the minimum processed B
                        Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end)); % logic matrix
                        locmin_1_min = find(Bmin_MinPos~=0,1); % returns the null's index in Bmin_MinPos
                        locmin_1_min = round(length(Bmin)/2)+locmin_1_min-1; % returns the null's index in Bmin 

                        peak_compare((p+1),1) = p; % saves the extension factor
                        peak_compare((p+1),2) = max(Bmin(locmin_1_min:end)); % saves the max value which is the PSL height
                        peak_compare((p+1),3) = N;
                    end % ends ext factor (p) loop  

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                    % this chunk stores configs for the different p values (with the same M and L)
                    no_ofindices = find(peak_compare(:,2)<=-13); % saves indices of successful configs
                    Configs_good(1:length(no_ofindices),1) = M; 
                    Configs_good(1:length(no_ofindices),2) = peak_compare(no_ofindices,1);
                    Configs_good(1:length(no_ofindices),3) = M*(peak_compare(no_ofindices,3)-1)+1;
                    Configs_good(1:length(no_ofindices),4) = M*((peak_compare(no_ofindices,1))+1); % stores totalM length

                    % concatenates the cases stored above to the previously tested success cases
                    finalConfigs_good = cat(1,finalConfigs_good(:,:), Configs_good(:,:));

                    no_ofindices_bad = find(peak_compare(:,2)>-13);
                    Configs_bad(1:length(no_ofindices_bad),1) = M;
                    Configs_bad(1:length(no_ofindices_bad),2) = peak_compare(no_ofindices_bad,1);
                    Configs_bad(1:length(no_ofindices_bad),3) = M*(peak_compare(no_ofindices_bad,3)-1)+1;
                    Configs_bad(1:length(no_ofindices_bad),4) = M*((peak_compare(no_ofindices_bad,1))+1);

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
        
        % analyzes configurations and selects best option
        function SelectConfiguration(flag, sensors) 
        % if flag == 0 simulated data, flag == 1 real data
        % sensors is an array of 1s and 0s to simulate on and off sensors

            close all;
            load('sub2s.mat'); 
            load('finalConfigs_good.mat');
            load('sub2s_bad.mat'); 
            load('finalConfigs_bad.mat');

            ElementsInUse = zeros(1,length(sensors));
            f = find(diff([0,sensors,0]==1)); % index of every change
            start1s = f(1:2:end-1);  % Start indices of 1s
            count1s = f(2:2:end)-start1s;  % Consecutive 1s counts
            [maxM,position] = max(count1s); %returns [max consecutive 1s, index in x]
            maxMstart = start1s(position);
            if maxM <2
                error('M is too small. This array cannot function as an NSA.')
            end

            first = start1s(1); 
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
                                for j = 1:(length(sensors)-length(sub2s_check)+1)
                                    if sum(lt(sensors((j):(j+length(sub2s_check)-1)),sub2s_check)) == 0
                                        match = 1
                                        sub2start = j
                                        M = finalConfigs_good(i,1);
                                        p = finalConfigs_good(i,2);
                                        ConfigPick = finalConfigs_good(i,:);
                                        %l = finalConfigs_good(i,3); %aperture

                                        break % out of j = 0:(length(sensors)-length(sub2s_check)) 
                                    end 
                                end
                                if match == 1
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
                for L = L_orig:-1:3  % check aperture decreasing
                    disp('Im checking a new aperture, bad'); disp(num2str(L))
                    sub2aperL = find(cell2mat(sub2s_bad(:,2))==L);  % saves indices of sub2s with aperture L
                    if isempty(sub2aperL) == 0 % there is a valid config of that L
                        lookhere = find(finalConfigs_bad(sub2aperL(1):sub2aperL(end),4)<=maxM) + sub2aperL(1)-1; % with valid sub1+ext
                        if isempty(lookhere)==0  % there is a valid config of that L (lookhere is not empty) 
                            for i = lookhere(1):lookhere(end)  % search through the valid sub2s 
                                sub2s_check = sub2s_bad{i,1};  % saves current sub2 being checked 
                                for j = 1:(length(sensors)-length(sub2s_check)+1)
                                    if sum(lt(sensors((j):(j+length(sub2s_check)-1)),sub2s_check)) > 0
                                        badcase = 1;
                                        sub2start = j;
                                        M = finalConfigs_bad(i,1);
                                        p = finalConfigs_bad(i,2);
                                        ConfigPick = finalConfigs_bad(i,:);
                                        %l = finalConfigs_bad(i,3); %aperture

                                        break % j = 0:(length(sensors)-length(sub2s_check))
                                    end
                                    if badcase == 1
                                        break % out of i = lookhere(1):lookhere(end)
                                    end
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

            %%%% Plot and selecting sensors!!
            if (match == 1)||(badcase == 1)
                ElementsInUse(maxMstart:(maxMstart+ConfigPick(1,4)-1)) = 1;
                ElementsInUse(sub2start:ConfigPick(1,1):(sub2start+length(sub2s_check)-1)) = 1;
                if flag == 1   % real data
                    % ArraySystem.TemporalSpatial(data,filters,processor,plots,M,p,L)
                    ArraySystem.TemporalSpatial(0,1,2,1,M,p,L); %figure;
                    % ArraySystem.NestedBP(M, p, L, plots)
                    ArraySystem.NestedBP(M,p,L,1);
                elseif flag ==0  % simulated data
                    % ArraySystem.TemporalSpatial(data,filters,processor,plots,M,p,L)
                    ArraySystem.TemporalSpatial(1,0,2,1,M,p,L); %figure;
                    % ArraySystem.NestedBP(M, p, L, plots)
                    ArraySystem.NestedBP(M,p,L,1);
                end
            elseif badcase == 1
                msgbox('This configuration is NOT guaranteed to accurately and adequately detect the signal.','Warning','warn')
            else
                msgbox('Array cannot be analyzed as an NSA.', 'Error', 'error')
            end
        end
    end
end
