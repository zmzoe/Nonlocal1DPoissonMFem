


%% ureal

ureal={@(x) (-x+1/4);
       @(x) (x-3/4)};
% (0,1/2)
ureala=@(x) -x+1/4;
% (1/2,1)
urealb=@(x) x-3/4;


preal={@(x) (-1);
       @(x) (1)};
%% preal 
% (0,1/2)
preala=@(x) -1;
% (1/2,1)
prealb=@(x) 1;


%% f
if delta>1/4
f1={@(x) ((-8/(3*delta^4))*x^3 + (4/delta^4)*x^2 + (-(2*(3/(2*delta^2) - 4))/delta^2)*x + (2*(delta - (delta/2 + (delta*(2*delta - 1))/2)/delta + 1/(3*delta^2) - 2))/delta^2);
    @(x) (-4/(3*delta^4))*x^3 + (4/delta^4)*x^2 + (-(2*(3/(2*delta^2) - 2))/delta^2)*x + (2*((4*delta)/3 + 1/(3*delta^2) - 2))/delta^2;
    @(x) (8/(3*delta^4))*x^3 + (-(2*((2*delta + 2)/delta^2 - 2/delta + 2/delta^2))/delta^2)*x^2 + ((2*((4*delta + 3/2)/delta^2 - (4*delta + 4)/delta + 2/delta^2))/delta^2)*x - (2*(delta + ((3*delta)/2 + 1/3)/delta^2 - ((3*delta)/2 + ((2*delta + 1)*(delta + 2))/2 + 1/2)/delta + 2/(3*delta^2)))/delta^2;
    @(x) (4/(3*delta^4))*x^3 + (-4/delta^4)*x^2 + ((2*(2/delta^2 - 2))/delta^2)*x - (2*((4*delta)/3 + 2/(3*delta^2) - 2))/delta^2};
else
f2={@(x) (-4/(3*delta^4))*x^3 + (4/delta^2)*x - (2*(delta/3 + (delta/2 + (delta*(2*delta - 1))/2)/delta))/delta^2;
    @(x) (-4/(3*delta^4))*x^3 + (4/delta^4)*x^2 + (-(2*(3/(2*delta^2) - 2))/delta^2)*x + (2*((4*delta)/3 + 1/(3*delta^2) - 2))/delta^2;
    @(x) (4/(3*delta^4))*x^3 + (-(2*((2*delta + 2)/delta^2 - 2/delta))/delta^2)*x^2 + ((2*((4*delta + 3/2)/delta^2 - (4*delta + 4)/delta + 2))/delta^2)*x + (2*(delta/3 - ((3*delta)/2 + 1/3)/delta^2 + ((3*delta)/2 + ((2*delta + 1)*(delta + 2))/2 + 1/2)/delta - 2))/delta^2;
    @(x) (4/(3*delta^4))*x^3 + (-4/delta^4)*x^2 + ((2*(2/delta^2 - 2))/delta^2)*x - (2*((4*delta)/3 + 2/(3*delta^2) - 2))/delta^2};

end

%(0,delta)
F1a=@(x) (-8/(3*delta^4))*x^3 + (4/delta^4)*x^2 + (-(2*(3/(2*delta^2) - 4))/delta^2)*x + (2*(delta - (delta/2 + (delta*(2*delta - 1))/2)/delta + 1/(3*delta^2) - 2))/delta^2;
F1b=@(x) (-4/(3*delta^4))*x^3 + (4/delta^2)*x - (2*(delta/3 + (delta/2 + (delta*(2*delta - 1))/2)/delta))/delta^2;

%(delta,1/2)
F2=@(x) (-4/(3*delta^4))*x^3 + (4/delta^4)*x^2 + (-(2*(3/(2*delta^2) - 2))/delta^2)*x + (2*((4*delta)/3 + 1/(3*delta^2) - 2))/delta^2;

%(1/2,delta+1/2)
F3a=@(x) (8/(3*delta^4))*x^3 + (-(2*((2*delta + 2)/delta^2 - 2/delta + 2/delta^2))/delta^2)*x^2 + ((2*((4*delta + 3/2)/delta^2 - (4*delta + 4)/delta + 2/delta^2))/delta^2)*x - (2*(delta + ((3*delta)/2 + 1/3)/delta^2 - ((3*delta)/2 + ((2*delta + 1)*(delta + 2))/2 + 1/2)/delta + 2/(3*delta^2)))/delta^2;
F3b=@(x) (4/(3*delta^4))*x^3 + (-(2*((2*delta + 2)/delta^2 - 2/delta))/delta^2)*x^2 + ((2*((4*delta + 3/2)/delta^2 - (4*delta + 4)/delta + 2))/delta^2)*x + (2*(delta/3 - ((3*delta)/2 + 1/3)/delta^2 + ((3*delta)/2 + ((2*delta + 1)*(delta + 2))/2 + 1/2)/delta - 2))/delta^2;

%(delta+1/2,1)
F4=@(x) (4/(3*delta^4))*x^3 + (-4/delta^4)*x^2 + ((2*(2/delta^2 - 2))/delta^2)*x - (2*((4*delta)/3 + 2/(3*delta^2) - 2))/delta^2;



%%




g={@(x) (-1-((-2/delta^2)*x^2 + (4/delta)*x - (2*(delta/4 + (delta*(2*delta - 1))/4))/delta^2));
   @(x) (-2);
   @(x) 1-((2/delta^2)*x^2 + (-(2*(2*delta + 2))/delta^2)*x + (2*((3*delta)/4 + ((2*delta + 1)*(delta + 2))/4 + 1/4))/delta^2);
   @(x) 2};

%(0,delta)
g1=@(x)  -1-((-2/delta^2)*x^2 + (4/delta)*x - (2*(delta/4 + (delta*(2*delta - 1))/4))/delta^2);

%(delta,1/2)
g2=@(x)  -2;

%(1/2,delta+1/2)
g3=@(x)  1-((2/delta^2)*x^2 + (-(2*(2*delta + 2))/delta^2)*x + (2*((3*delta)/4 + ((2*delta + 1)*(delta + 2))/4 + 1/4))/delta^2);

%(delta+1/2,1)
g4=@(x)  2;




%%


HH={@(x) 0;
    @(x) -2/(delta^2)*(2*delta+2*x-1);
    @(x) 0;
    @(x) -2/(delta^2)*(-2*delta+2*(1-x))};
%(0,delta)
HH1=@(x)  0;

%(delta,1/2)
HH2=@(x)  -2/(delta^2)*(2*delta+2*x-1);

%(1/2,delta+1/2)
HH3=@(x)  0;

%(delta+1/2,1)
HH4=@(x)  -2/(delta^2)*(-2*delta+2*(1-x));



