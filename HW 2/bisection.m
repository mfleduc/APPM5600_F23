function [root, cnt, errNum, errMsg] = bisection(fHandle, interval, tol, maxIters)
%BISECTION Implementation of bisection method
% Error messages:
%   0:Success
%   1:Root not found due to no detected sign change on interval
%   2:Root not found, max iterations reached
a = interval(1);b = interval(2);%
fa = fHandle(a);
fb = fHandle(b);
if fa*fb>0 %Handling the case where sign(a) and sign(b) are different
    newPts = linspace(a, b, 5);%checking if the initial choice of a and b 
                                %were just bad
    for ii = 1:length(newPts)-1
        if fHandle(newPts(ii))*fHandle(newPts(ii+1))<0
            a = newPts(ii);%If we find a sign change in one of the new 
                            %intervals we are cooking
            b = newPts(ii+1);
            fa = fHandle(a);
            fb = fHandle(b);
            rootExists = 1;
            break
        else
            rootExists = 0;
        end
    end
    if ~rootExists
        %If we don't find a sign change, output the midpoint
        root = (a+b)/2;
        cnt = 0;
        errNum = 1;
        errMsg = "Root not found: No sign change detected";
        return
    end
end
del = tol+1;
cnt = 0;
while del>tol && cnt<maxIters%Bisection method 
    cnt = cnt+1;
    c = (a+b)/2;
    fc = fHandle(c);
    if fa*fc<0
        b = c;
        fb = fc;
    else
        a = c;
        fa = fc;
    end
    del = abs( (a-b)/b );
end
root = (a+b)/2;
if cnt<maxIters
    errNum = 0;
    errMsg = "Successful";
else
    errNum = 2;
    errMsg = "Root not found: Maximum iterations reached";
end
end

