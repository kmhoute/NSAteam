classdef BP_Formation
    methods(Static)
        % creates the Nested Sensor Array Beampattern
        function [u, B1, B2, Bmin, Bprod, N, add1] = Nested(M, p, L, plots)
        %%% Generates the nested beam pattern
        %%% plots = 1, a plot is created
        %%% plots = 0, a plot is not created
        
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
            %%%%The first argument in BeampatternLinearArray is set to 0 so that
            %%%%that function doesn't plot graphs. The second parameter is the
            %%%%number of sensors in a subarray. The third parameter is the number
            %%%%of extra sensors in a subarray for extension. The fourth parameter
            %%%%is the intersensor spacing in terms of lambda. So, need to specify
            %%%%0.5 as well.
            B1 = BP_Formation.Linear(0, M+add1, 0, U1*0.5);
            B1 = B1/max(abs(B1));
            B2 = BP_Formation.Linear(0, N+add2, 0, U2*0.5);
            B2 = B2/max(abs(B2));


            Bmin = min(abs([B1;B2]));
            Bmin = 20*log10(abs(Bmin));

            Bprod = B1.*conj(B2);
            Bprod = 10*log10(abs(Bprod));
            
            
            if plots == 1
                f1 = figure;
                a1 = axes('Parent', f1, 'FontSize', 16, 'FontWeight', 'Bold');%...
                    %'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
                hold(a1, 'all');
                box(a1, 'on');
                grid(a1, 'on');
                plot(a1, u, 20*log10(abs(B1)), 'LineStyle', '-.', 'LineWidth', 1, 'Color', 'Blue');
                hold(a1, 'all');
                plot(a1, u, 20*log10(abs(B2)), 'LineStyle', '-.', 'LineWidth', 3, 'Color', 'Red');
                hold(a1, 'all');

                plot(a1, u, Bmin, 'LineWidth', 2, 'Color', 'Black');
                hold on;
                plot(a1, u, Bprod, '--','LineWidth', 1, 'Color', [0 0.6 0]);
                hold on;
                plot(a1, [-1 1],[-13 -13],'--','LineWidth',1,'Color','k');
                xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
                ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
                % set(gcf, 'Position', [50 1 1317 689]); % Maximicoprimepairs.pdfze figure.         
                legend('Subarray 1','Subarray 2', 'Min', 'Product');
                xlim([-1 1]);
                ylim([-30 0]);
                title(['add1: ',num2str(add1),' and U2: ',num2str(U2)])
            elseif (plots ~= 0)|(plots ~= 1)
                disp('Error; plots must be 0 or 1')
            end
        end
        
        % creates the Uniform Linear Array Beampattern
        function [B, u] = Linear(flag, N, varargin)
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
            for idx = 1:length(u);
                v = exp(1i*(n-(N-1)/2)*2*pi*(u(idx)-shift)*d);
            %   v = exp(1i*(n)*2*pi*(u(idx)-shift)*d);
                B(idx) = conj(w)*v;
            end
            % Normalize by the absolute value at u = 0
            B = B/sum(w);
            if flag
                f = figure;
                set(gcf, 'Position', [50 1 1317 689]); % Maximize figure.     
                a = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16, ...
                    'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
                plot(a, u, (real(B)), 'LineWidth', 3, 'Color', 'Black');
                hold on;
                plot(a, u, (imag(B)), 'LineWidth', 3, 'Color', 'red');
            %   title('Beam Pattern in U - Space', 'FontWeight', 'Bold', ...
            %         'FontSize', 16);
                grid on;
                xlabel('u', 'FontSize', 16, 'FontWeight', 'Bold');
                ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
            %   ylim([-30 0]);
                xlim([-1 1]);
        end
end
    end
end
