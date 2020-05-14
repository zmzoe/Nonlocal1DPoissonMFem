

% Input : Alpha  the transformation matrix
% Output : The matrix of G_delta^* Psi 
%          
%
% II : G psi over II-th interval
%
% aa_i x^2 + bb_i x + cc_i
%---------------------------------------%

function [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,II)


                [alpha_i,beta_i,gamma_i]=gradp_dg(delta,m,II);
                
                 aa_i=alpha_i*Alpha';
                 bb_i=beta_i*Alpha';
                 cc_i=gamma_i*Alpha';

end