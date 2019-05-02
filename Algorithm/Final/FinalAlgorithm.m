%Final Algorithm
function FinalAlgorithm(real_or_simulated,varargin)
tic
% real_or_simulated 
    % 1 = real data
    % 0 = simulated data
        % Simulated data must provide 
        % EITHER values of M, p, L
        % OR one matrix of 1/0 in brackets %ex. [1 1 1 0 0 1]
%
        
close all 
clc
if real_or_simulated == 1 %use real data
% Actual Data
    [filename, path] = uigetfile();                                            % Get the name and the path of the file 
    load([path filename]);                                                     % Load the file.
    totalData = Rearrange_Data(totalData);
    save('totalData.mat')
    % PSD analysis 
        shape_of_window_kaiser = 0.05;                                         % Shape of the kaiser window  
        sample_rate = 44100;
        w = kaiser(length(totalData),shape_of_window_kaiser);                  % Kaiser windowing
        rbw = enbw(w,sample_rate);                                             % Returns the two-sided equivalent noise bandwidth, bw, in Hz.
            for i = 1:size(totalData,2)                                        % for loop to calculate each microphone PSD by it self   
                [ sxx(:,i),F(:,i) ] = periodogram(totalData(:,i)-mean(totalData(:,i)),w,length(totalData),sample_rate,'psd');   
                                                                               % sxx = the PSD and F = frequencies specified in the vector.
                SNR_1(:,i) = snr(sxx(:,i),F(:,i),rbw,'power');                 % Singal to nouse ration for each Sensor   
            end
        sxxx = mean(sxx,2);                                                     % the avrage PSD for all of the sensor  
        ff = mean(F,2);                                                        % the avaege frequencies specified for all of the sensor       
        SNR_check = snr(sxxx,ff(:,1),rbw,'power');                             % singal to noise ratio and the PSD OF the signal
        %snr(sxxx,ff,rbw,'power');                                               % Display the figure
        %set(gca,'FontSize',15);
        threshold_SNR= mean(SNR_1,2,'omitnan');
        var_1 = var(SNR_1,'omitnan');                                                    % Variance of the signal to noise values 
            if var_1 < 0.5
               sensors = ones(1,size(totalData,2));  
            else 
               threshold = threshold_SNR - sqrt(var_1);   
                 for j = 1:size(totalData,2)                                   % For loop to check if the resulten SNR is higher than the SNR...
                                                                               % of the individule sensor
                     if SNR_1(:,j) > threshold                                 % If it is higher 
                        sensor_on(:,j) = 1;                                    % Then place 1 
                     else                                                      % If it is not 
                        sensor_on(:,j) = 0;                                    % Then place 0
                     end
                 end   
            end
        good_sensor_data = totalData.*sensor_on;                               % actual good data
        sensors = sensor_on(1,:); 
        figure
        plot(SNR_1,'o','lineWidth',2)
        set(gca,'FontSize',15)
        hold on 
        yline(threshold_SNR,'r','lineWidth',2)
        yline(threshold,'k','lineWidth',2)
        xlabel('\bf Sensors','FontSize',20,'FontWeight','bold','Color','k');
        ylabel(' \bf Signal to noise ratio dB ','FontSize',20,'FontWeight','bold','Color','k');
        title(' \bf Determination of valid sensors ','FontSize',20,'FontWeight','bold','Color','K');
    %     annotation('textarrow','String',' Average of SNR','FontSize',20);
    %     annotation('textarrow','String',' Lower bound \sigma^2','FontSize',10);
        legend('Sensor','Average of SNR','Lower bound \sigma^2','FontSize',20,'Location','southeast');

        MainAlgorithm(1,sensors);
        
elseif real_or_simulated == 0 %use simulated data
%     if length(varargin) == 2
%         error('Incorrect inputs')
%         disp('Varargin must have 1 or 3 elements')
%     end
    vars = cell2mat(varargin);
    if length(varargin) == 3
    vars = cell2mat(varargin);    
    M = vars(1); 
    p = vars(2);
    L = vars(3);
    NSA_Analysis(1,0,2,1,M,p,L);
    BP_Formation.Nested(M,p,L,1);
    elseif length(varargin) == 1
        vars = cell2mat(varargin);
        MainAlgorithm(0,vars);
    end
end
toc
end
