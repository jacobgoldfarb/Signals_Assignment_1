
convolute([1,3,2,1,2,2,5],[1,0,1,0,1])
test();
function test()
    for i = 0:1:1000
       x1 = fix(rand([1 fix(100 * rand)]) * 100);
       x2 = [1,10000,1,0,1];% fix(rand([1 fix(100 * rand)]) * 100);
       
       target = conv(x1, x2);
       actual = cast(convolute(x1, x2), 'double');
       delta = actual - target
       assert(all(delta == 0));
    end
end

function conv = convolute(x, h)
    
    m = length(x);
    n = length(h);
    % pad x and h with zeros
    x = [x, zeros(1,n)];
    h = [h, zeros(1,m)];
    
    % for each element in x, multiply it by every element in h and append
    % the sum to y.
    for i = 1:(m + n - 1)
        sum = 0;
           for j = 1:m
               if (i - j + 1 > 0)
                sum = sum + h(i - j + 1) * x(j);
               end
           end
        y(i) = sum;
    end
    conv = y;
end
