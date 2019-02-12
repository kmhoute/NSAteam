N=64;

load 'Senior Design\Measurement_Data\3\90'
    [totalData order]=Rearrange_Data(totalData);
    %Apply a highpass filter to remove noise
    band = designfilt('highpassiir','FilterOrder',50,...
    'PassbandFrequency', 3500, 'PassbandRipple',.2,...
    'SampleRate',44100); %This filter passes frequencies above 4000
    totalData=filter(band,totalData);
    
    
    totalData=totalData(100:44200,3:N);
time = linspace(0,1,length(totalData));    
%subplot(2,1,1)
plot(time,totalData(:,1));
%ylim([-.005,.005]);
%xlim([.015,.02])
hold on
%subplot(2,1,2)
plot(time,totalData(:,15));
plot(time,totalData(:,29));
plot(time,totalData(:,33));
plot(time,totalData(:,47));
plot(time,totalData(:,61));

ylim([-.002,.002]);
xlim([.018,.02])
legend