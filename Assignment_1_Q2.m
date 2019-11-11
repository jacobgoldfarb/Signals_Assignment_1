 
sample_rate = 0.01
t = 0.0:sample_rate:15;
xt = exp(-2*pi*j*t)%1+cos(2*pi*t)/4 + cos(2*pi*t*2)/2 + cos(2*pi*t*3)/3;


[ak, k] = fourier_series(xt, t);
[retrieved_xt, retreived_t] = inv_fourier_series(ak, k);


function [dv, iv] = fourier_series(xt, t)
   figure('Name', 'Inputted x(t)'), plot(t,xt);
   figure('Name', 'Inputted x(t) (imaginary part)'), plot(t,imag(xt));
   
   %dt is the difference between two x values.
   dt = t(2) - t(1);
   %the period, and number of samples in the period (N) is retrieved.
   [period, N] = get_ft(xt, t);
   omega = (2 * pi) / period;
   
   %The signal is integrated from dt to period, stepping by dt each time.
   integration_interval = dt : dt : period;
   
   % k is an arbitrary range, in this case from -20 to 20.
   k_min = -20;
   k_max = 20;
   ks = k_min:1:k_max;
   ak = zeros(1, length(ks)); % instantiate array of aks as zeroes
   
   % Only x(t) values occuring in the fundamental period are integrated
   % over.
   xt_range = xt(1:N)
   figure('Name', 'Fundamental Period'), plot(integration_interval, xt_range)
   
   % a(k) values are assigned iteratively for each k in the arbitrary
   % range.
   for k = ks
       % element-wise vector multiplication and vector summation is used to determine the
       % integral for a single k value.
       integral = sum((xt(1: N) .* exp(-1i * omega * k * integration_interval ) * dt));
       % the integral is then divided by the period and assigned to a(k)
       ak(k + k_max + 1) = integral / period;
   end
   % return values are assigned 
   dv = ak;
   iv = ks;
   figure('Name', 'ak real'), plot(ks,real(ak));
   figure('Name', 'ak imaginary'), plot(ks,imag(ak));
end

function [dv, iv] = inv_fourier_series(ak, ks)
    %some arbitrary t array is allocated.
    ts = 0:0.01:5;
    %omega is hardcoded
    omega = 2 * pi;
    xts = []; % instantiate empty array of xts
    for singular_t = ts
        % summation occurs for every t in x(t), similar to the fourier
        % series function.
        xts_sum_for_one_t = sum(ak .* exp(-1i * omega * singular_t * ks));
        xts = [xts xts_sum_for_one_t];  
    end
    %return variables are assigned
    dv = xts;
    iv = ts;
    figure('Name', 'retreived X(t) vs t'), plot(ts,xts);
end

%The signal must contain at least two periods
function [period, period_index_count] = get_ft(xt, t)
    %Periods will only be considered if they consist of at least 3 samples
    MIN_PERIOD_SAMPLES = 3;
    %Some tolerance is used to compare sections of the signal
    TOLERANCE = 0.15;
    %The entire left half of the signal is compared to the right half of the signal    
    left_period = xt(1 : length(xt)/2);
    right_period = xt(length(xt)/2 : length(xt)/2 + length(left_period) - 1);
    %default values
    period = -1; 
    period_index_count = 0;
    uniques_in_period = 0;

    while length(left_period) > 1%MIN_PERIOD_SAMPLES
        %In every iteration, the righter most element of the left period is
        %moved to the right half, and the right most element of the right
        %period is removed.
        right_element_of_left_period = left_period(end);
        left_period(end) = [];
        right_period = [right_element_of_left_period right_period];
        right_period = right_period(1: length(right_period) - 2);
        %The resulting left half and right half are compared for equality
        %(against some tolerance).

        if (all(abs(left_period - right_period) < TOLERANCE) == 1)
            %If the signal sections are the same, then one section must 
            %also contain each unique value from the last calculated
            %period.
            if period == -1 || length(uniquetol(real(left_period))) + length(uniquetol(imag(left_period))) ==  uniques_in_period
                uniques_in_period = length(uniquetol(real(left_period))) + length(uniquetol(imag(left_period)));
                period = t(length(left_period) + 1) - t(1);
                period_index_count = length(left_period);
            end
        end
    end
end
