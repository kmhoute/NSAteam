% psl
psl = zeros(1, 15)
success = zeros(1, 15)
M_length = [1:15]
for m = 1:15
    psl(m) = 20*log10(abs([sin(0.5*pi*3)]/[m*sin(0.5*pi*3/m)]));
    if psl(m) <= -13
        success(m) = psl(m);
    else
        success(m) = 0;
    end
end

psl
success
% plot(M_length,success, 'bx'); hold on;
% plot(M_length,psl, 'rx')

M = 12
N = 5
u = -1:.0003/(M*N):1;
sub1 = 20*log10(abs([sin(0.5*pi*M.*u)]./[M*sin(0.5*pi.*u)]));
sub2 = 20*log10(abs([sin(0.5*pi*M*N.*u)]./[N*sin(0.5*pi*M.*u)]));
plot(u,sub1, 'b-');hold on;
plot(u,sub2, 'r--');
yline(-13)
xline(3/(M*N))
xline(3/((25/24)*M*N))
ylim([-35 0])

%% fudge factor %
data = [8, 9, 0.0399, 0.9576;...
        8, 5, 0.0725525, 0.96736667;...
        4, 9, 0.0798, 0.9576;...
        5, 12, 0.04779, 0.9558;...
        12, 5, 0.04837, 0.9674]
                  
plot(data(:,2),data(:,4), 'bo'), hold on;
% plot(data(:,2),data(:,1), 'rx')
ylim([0.94 1])
