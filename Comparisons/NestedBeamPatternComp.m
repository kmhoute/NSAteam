thresh = 4; % L is the min # of sensors to be tested in the for loop and is 
            % therefore the length of subarray 1. This means that the first
            % iteration of the for loop is full array analysis

M = thresh; % + add1=p*M
p = 0;
peak_compare = zeros(64, 4);
for L = (thresh*3):64  % comment this line to run beampattern only; 
                  % uncomment for prod vs min analysis on # of sensors
    [u, B1, B2, Bmin, Bprod, N, add1] = ProdMinBPComp(M,p,L);
    
    
     Bmin_MinPos = islocalmin(Bmin(floor((length(Bmin)/2)):end));
     Bprod_MinPos = islocalmin(Bprod(floor((length(Bprod)/2)):end));
    % [peaks_B1, ~]    = findpeaks(B1, u,'SortStr', 'descend');
    % [peaks_B2,~]    = findpeaks(B2, u, 'SortStr', 'descend');
    peaks_Bmin  = max(Bmin(Bmin_MinPos:end));
    peaks_Bprod = max(Bprod(Bprod_MinPos:end));
    
%     peak_compare(L, 1) = peaks_B1(2);
%     peak_compare(L, 2) = peaks_B2(1);
    peak_compare(L, 3) = peaks_Bmin;
    peak_compare(L, 4) = peaks_Bprod;
end  %%%%% comment this line and the following lines to run beampattern only; 
       %%%% uncomment for prod vs min analysis on # of sensors

sensors = (1:64)';
% sub1 = peak_compare(:,1);
% sub2 = peak_compare(:,2);
prod = peak_compare(:,3);
min = peak_compare(:,4);
% scatter(sensors, sub1, 'b+')
% hold on;
% scatter(sensors, sub2, 'rx')
% hold on;
scatter(sensors, prod, 'k')
hold on;
scatter(sensors, min, 'g')
legend('product', 'minimum') % 'subarray 1', 'subarray 2', 
ylabel('peaks','FontWeight','bold');
xlabel('sensors used','FontWeight','bold');
title('subarray 1 length = ', num2str(M), ' and subarray 2 length = ',num2str(N), 'FontWeight','bold');

