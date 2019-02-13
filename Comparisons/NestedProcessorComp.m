load 't75.mat'
L = 16; % L is the min # of sensors to be tested in the for loop and is 
        % therefore the length of subarray 1. This means that the first
        % iteration of the for loop is full array analysis
 % [totalData, order]=Rearrange_Data(totalData); % only use when running
                                                % gathered data
%%% Apply a highpass filter to remove noise
band = designfilt('highpassiir','FilterOrder',50,...
'PassbandFrequency', 3500, 'PassbandRipple',.2,...
'SampleRate',44100); %This filter passes frequencies above 4000
totalData=filter(band,totalData);

peak_compare = zeros(64, 2);
% totalData=totalData(100:44200,1:64); % only use when running
                                       % gathered data
N = 64;
for N = (L):64  % comment this line to run beampattern only; 
                  % uncomment for prod vs min analysis on # of sensors
    sparseindices1 = (1:L); %%%We add 1 at the end because MATLAB...
                                   %%%index starts at 1, not 0
    sparseindices2 = (1:4:(N)); %%%We add 1 at the end because MATLAB...
                                   %%%index starts at 1, not 0
    data1 = zeros(size(totalData));
    data2 = zeros(size(totalData));
    data1(:,sparseindices1) = totalData(:,sparseindices1);
    data2(:,sparseindices2) = totalData(:,sparseindices2);

    for processor = 1:2
        if (processor==1)
            F = fft(data1,64*10,2) .* conj(fft(data2,64*10,2)); %Find the product
            mF = mean(abs(F)); %Find the average
            mF = 10*log10(fliplr(fftshift(mF)/max(abs(mF)))); %Normalize, fftshift, flip L to R & convert to dB
        elseif (processor==2)
            F = min(abs(fft(data1,64*10,2)) , abs(fft(data2,64*10,2))); %Find the minimum of the absolute values
            mF = mean(abs(F)); %Find the average
            mF = 20*log10(fliplr(fftshift(mF)/max(abs(mF)))); %Normalize, fftshift, flip L to R & convert to dB
        else 
            disp('Error: Input argument must be 1 or 2.');
            return;
        end
        w = linspace(-1,1,length(mF));
        [peaks, loc] = findpeaks(mF, w,'SortStr', 'descend');
        
        if (processor==1)
            peak_compare(N, 1) = (peaks(3)+peaks(4))/2;
        else
            peak_compare(N, 2) = (peaks(3)+peaks(4))/2;
        end

        %%%%%%%%%%%Plot the spatial spectrum
     %%%% uncomment this line and the following lines to run beampattern only; 
     %%%% comment for prod vs min analysis on # of sensors

%         if (processor==1)
%             plot(acosd(w),mF, 'r','LineWidth',1.5);
%         else
%             plot(acosd(w),mF, 'b-.','LineWidth',1);
%         end
%         ylabel('Power dB','FontWeight','bold');
%         xlabel('cos(\theta)','FontWeight','bold');
%         title('Spatial Spectral Estimation: 75 degrees, subarray lengths 16 and 4','FontWeight','bold');
%         xlim([0 180]);
%         legend('product', 'min');
%         ylim([-20 0]);

        hold on;
    end
end  %%%%% comment this line and the following lines to run beampattern only; 
       %%%% uncomment for prod vs min analysis on # of sensors

sensors = (1:64)';
product = peak_compare(:,1);
minimum = peak_compare(:,2);
scatter(sensors, product, 'b')
hold on;
scatter(sensors, minimum, 'r')
legend('product', 'minimum')
ylabel('peaks','FontWeight','bold');
xlabel('sensors used','FontWeight','bold');
title('75 degrees, subarray lengths 16 and 2','FontWeight','bold');

