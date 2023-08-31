%% HW 1 Problem 2
close all; clear variables

a0 = 4.8;
b0 = 5.31;
tol = 1E-4;
%% Part 1: Factored 
del = 1;
count = 0;
a=a0;b=b0;
fa = (a0-5)^9;
fb = (b0-5)^9;
while del>tol && count<1e5
    c = ( a+b )/2;
    fc = (c-5)^9;
    switch sign( fc )
        case sign(fa)
            a = c;
            fa = ( a-5 )^9;
        case sign(fb)
            b = c;
            fb = ( b-5 )^9;
        otherwise
            keyboard
    end
    del = abs( b-a );
    count = count+1;
end
root = (a+b)/2;
fprintf( ['Expanded form: Root found at apprroximately x=%.5f\n',...
 'Iterations=%d\n\n'], root,count );
%% Part 2
del = 1;
count = 0;
a=a0;b=b0;
fa = f1(a0);
fb = f1(b0);
while del>tol && count<1e5
    c = ( a+b )/2;
    fc = f1(c);
    switch sign( fc )
        case sign(fa)
            a = c;
            fa = f1(c);
        case sign(fb)
            b = c;
            fb = f1(c);
        otherwise
            keyboard
    end
    del = abs( b-a );
    count = count+1;
end
root = (a+b)/2;
fprintf( ['Expanded form: Root found at apprroximately x=%.5f\n',...
    'Iterations=%d\n\n'], root,count );

function y = f1(x)
y = zeros(size(x));
for ii = 0:9
    y = y+nchoosek( 9, ii )*(x).^ii*(-5)^( 9-ii );
end
end