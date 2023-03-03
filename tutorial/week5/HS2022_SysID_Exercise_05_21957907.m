function [R_u,lags,omegas,U_cos,U_cos_Wf,U_cos_Wt] = HS2022_SysID_Exercise_05_21957907()

%% Output format specification
% R_u must be a 3xN matrix
% lags must be a 1xN vector
% omegas must be a 1xN vector
% U_cos must be a 1xN vector
% U_cos_Wf must be a 1xN vector
% U_cos_Wt must be a 1xN vector
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[u_prbs,u_rand,u_cos] = HS2022_SysID_Exercise_05_GenerateData(LegiNumber);

N = length(u_prbs);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the input signals u_prbs, u_randn and u_cos to solve the problem.

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% 1. Calculation of autocorrelation

R_u = NaN * ones(3,N);

lags = NaN * ones(1,N);


%% 2. Smoothing

omegas = NaN * ones(1,N);

U_cos = NaN * ones(1,N);

U_cos_Wf = NaN * ones(1,N);

U_cos_Wt = NaN * ones(1,N);


%% Functions

    function [lags_w,wHann] = WtHann(gamma,N)
        %------------------------------------------------------------------
        %
        %   [lags,wHann] = WtHann(gamma,N)
        %
        %   Create a Hann window with width parameter gamma and data length N.
        %   The Hann window is a raised cosine.
        %
        %   Roy Smith,  18 October, 2017.
        %
        %------------------------------------------------------------------

        if nargin == 0,
            disp('Syntax: [lags,w] = WtHann(gamma,N)')
            return
        elseif nargin ~= 2,
            error('incorrect number of input arguments (2 expected)')
            return
        end

        %   basic parameter checking
        if length(gamma) > 1,
            error('Width parameter, gamma, must be a scalar');
        end
        if round(gamma) ~= gamma,
            error('Width parameter, gamma, must be an integer');
        end
        if gamma < 1,
            error('Width parameter, gamma, must be positive');
        end
        if length(N) > 1,
            error('Calculation length, N, must be a scalar');
        end
        if round(N) ~= N,
            error('Calculation length, N, must be an integer');
        end
        if N < 1,
            error('Calculation length, N, must be positive');
        end

        lags_w = [floor(-N/2+1):floor(N/2)]';
        wHann = 0*lags_w;
        idx = find(abs(lags_w) <= gamma);
        wHann(idx) = 0.5*(1+cos(pi*lags_w(idx)/gamma));

    end

%--------

    function [omega,WHann] = WfHann(gamma,N)
        %------------------------------------------------------------------
        %
        %   [omega,WHann] = WfHann(gamma,N)
        %
        %   Create a frequency domain Hann window with width parameter gamma
        %   and data length N.  The Hann window is a raised cosine.
        %
        %   Roy Smith,  18 October, 2017.
        %
        %                6 November, 2017.  Fixed bug in N even indexing.
        %
        %------------------------------------------------------------------

        if nargin == 0,
            disp('Syntax: [omega,W] = WfHann(gamma,N)')
            return
        elseif nargin ~= 2,
            error('incorrect number of input arguments (2 expected)')
            return
        end

        %   basic parameter checking
        if length(gamma) > 1,
            error('Width parameter, gamma, must be a scalar');
        end
        if round(gamma) ~= gamma,
            error('Width parameter, gamma, must be an integer');
        end
        if gamma < 1,
            error('Width parameter, gamma, must be positive');
        end
        if length(N) > 1,
            error('Calculation length, N, must be a scalar');
        end
        if round(N) ~= N,
            error('Calculation length, N, must be an integer');
        end
        if N < 1,
            error('Calculation length, N, must be positive');
        end

        %   The simplest approach is to define the window in the time domain and
        %   then transform it to the frequency domain.

        lags_w = [floor(-N/2+1):floor(N/2)]';
        wHann = 0*lags_w;
        idx = find(abs(lags_w) <= gamma);
        wHann(idx) = 0.5*(1+cos(pi*lags_w(idx)/gamma));

        %
        zidx = find(lags_w==0);    % index of the zero point.

        wH_raw = fft([wHann(zidx:N);wHann(1:zidx-1)]);
        WHann(zidx:N) = wH_raw(1:N-zidx+1);  % shift +ve freq to end
        WHann(1:zidx-1) = wH_raw(N-zidx+2:N);% shift > pi freq to beginning
        WHann = real(WHann);   % should actually be real
        omega = 2*pi/N*lags_w;

    end

end
