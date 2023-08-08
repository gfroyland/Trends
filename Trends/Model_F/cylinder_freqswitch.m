function y=cylinder_freqswitch(x,delta,alpha1,alpha2)
%call with alpha1/alpha2 irrational approximately 5/2

y=[simoncircle(x(1),delta),mod(x(2)+alpha1*(0.5+0.5*tanh(40*(x(1)-.5)))+alpha2*(0.5-0.5*tanh(40*(x(1)-.5))),1)];

    
