function y=simoncircle(x,delta),

if x<0.25,
    y=2*x;
elseif (x>=0.25) & (x<0.75)
    y=mod(delta+2*(x-0.25),1);
else
    y=0.5+2*(x-0.75);
end

