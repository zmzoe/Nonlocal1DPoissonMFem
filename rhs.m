%RHS

%% ureal and preal piecewise and periodic
% (0,1/2)
% (1/2,1)
ureal={@(x) (-x+1/4);
       @(x) (x-3/4)};


preal={@(x) (-1);
       @(x) (1)};

%% Gu
% (0,delta)
% (delta,1/2)
% (1/2,1/2+delta)
% (1/2+delta,1)

Gu={@(x) (2*(delta*(x - 1/4) - (x*(2*x - 1))/4 + ((delta - x)*(2*x - 2*delta + 1))/4))/delta^2;
    @(x) 1;
    @(x) (2*(((2*x - 1)*(x - 1))/4 - delta*(x - 3/4) + ((delta - x)*(2*delta - 2*x + 1))/4))/delta^2;
    @(x) -1};

%% f
if delta>1/4

f={@(x) (2*(delta + (2*(delta*(x - 1/4) - (x*(2*x - 1))/4 + ((delta - x)*(2*x - 2*delta + 1))/4))/delta + (- (2*x^3)/3 + 2*delta*x^2)/delta^2 - ((2*x^3)/3 - x^2 + x/2 - 1/12)/delta^2 - 1))/delta^2;    
   @(x) (2*((4*delta)/3 + 2*x - ((2*x^3)/3 - x^2 + x/2 - 1/12)/delta^2 - 1))/delta^2;
   @(x) (2*((4*delta)/3 + 2*x - ((2*x^3)/3 - x^2 + x/2 - 1/12)/delta^2 - 1))/delta^2;
   @(x) (2*((2*(((2*x - 1)*(x - 1))/4 - delta*(x - 3/4) + ((delta - x)*(2*delta - 2*x + 1))/4))/delta - delta - (delta*(2*x^2 - 2*x + 1/2) - x/2 + x^2 - (2*x^3)/3 + 1/12)/delta^2 + ((2*x^3)/3 - 2*x^2 + 2*x - 2/3)/delta^2 + 1))/delta^2;
   @(x) -(2*((4*delta)/3 + 2*x - ((2*x^3)/3 - 2*x^2 + 2*x - 2/3)/delta^2 - 2))/delta^2;
   @(x) -(2*((4*delta)/3 + 2*x - ((2*x^3)/3 - 2*x^2 + 2*x - 2/3)/delta^2 - 2))/delta^2};

else 

f={@(x) -(2*(delta/3 + 2*x - (2*(delta*(x - 1/4) - (x*(2*x - 1))/4 + ((delta - x)*(2*x - 2*delta + 1))/4))/delta - (- (2*x^3)/3 + 2*delta*x^2)/delta^2))/delta^2;
   @(x) 0;
   @(x) (2*((4*delta)/3 + 2*x - ((2*x^3)/3 - x^2 + x/2 - 1/12)/delta^2 - 1))/delta^2;      
   @(x) (2*(delta/3 + 2*x + (2*(((2*x - 1)*(x - 1))/4 - delta*(x - 3/4) + ((delta - x)*(2*delta - 2*x + 1))/4))/delta - (delta*(2*x^2 - 2*x + 1/2) - x/2 + x^2 - (2*x^3)/3 + 1/12)/delta^2 - 1))/delta^2;
   @(x) 0;
   @(x) -(2*((4*delta)/3 + 2*x - ((2*x^3)/3 - 2*x^2 + 2*x - 2/3)/delta^2 - 2))/delta^2};
        
 
end
%% g=preal-Gu

g={ @(x) -1-((2*(delta*(x - 1/4) - (x*(2*x - 1))/4 + ((delta - x)*(2*x - 2*delta + 1))/4))/delta^2);
    @(x) -2;
    @(x) 1-((2*(((2*x - 1)*(x - 1))/4 - delta*(x - 3/4) + ((delta - x)*(2*delta - 2*x + 1))/4))/delta^2);
    @(x) 2};


%% HH=-G_delta^* preal
%(0,1/2-delta)
%(1/2-delta,1/2)
%(1/2,1-delta)
%(1-delta,1)

HH={@(x) 0;
    @(x) -2/(delta^2)*(2*delta+2*x-1);
    @(x) 0;
    @(x) -2/(delta^2)*(-2*delta+2*(1-x))};


