function [x , y , z]  = Metropolis_inactive_Diff( ...
...
x                , y            , z             , sign_siz    , ...
Number_particles , Delta_t      , b             , mu_prior_xy , ...
var_prior_xy     , mu_prior_z   , var_prior_z   , rc          , ...
zc               , D            )   
 
% %%%%%%%%% INPUT %%%%%%%%% %
     % observation                   :  Transformed signal through ANSCOMBE transform
     % mu_proposal                   :  Proposal parameter for MU           ( Gamma distribution)
     % mu_back_proposal              :  Proposal parameter for MU_BACK      ( Gamma distribution)
     % mu_alpha                      :  Alpha parameter of MU prior         ( Gamma distribution)
     % mu_beta                       :  Beta parameter of MU prior          ( Gamma distribution)
     % mu_back_alpha                 :  Alpha parameter of MU_BACK prior    ( Gamma distribution)
     % mu_back_beta                  :  Alpha parameter of MU_BACK prior    ( Gamma distribution)
     % sign_siz                      :  Length of the observation signal
     % Number_particles              :  Number of defined particles
     % mu_acceptance_rate            :  Previous acceptance rate 
     % signal                        :  Original signal
     % x                             :  Positions of the particles in two dimentions  x(t) = [(r(t)/wxy)^2;(z(t)/wz)^2]
     % MU_OLD                        :  Previous MU         ( Emission rate of the particles )
     % MU_BACK_OLD                   :  Previous MU_BACK        ( Background emission rate )
     % Delta_t                       :  Data acquesition time
     % b                             :  Loads for active particles
     
% %%%%%%%% OUTPUT %%%%%%%%% %
     % MU                            :  Updated MU 
     % MU_BACK                       :  Updated MU_BACK
     % mu_acceptance_rate            :  Updated acceptance rate
     
% Accept or reject the proposals
     for n=1:Number_particles
         if b(n)==0
            x(n,1)       = ((-1)^floor(2*rand()))*mu_prior_xy/rc + sqrt(var_prior_xy)*randn()   ;
            y(n,1)       = ((-1)^floor(2*rand()))*mu_prior_xy/rc + sqrt(var_prior_xy)*randn()   ;
            z(n,1)       = ((-1)^floor(2*rand()))*mu_prior_z/zc + sqrt(var_prior_z)*randn()   ;
            
%             x(3*n-1,1)       = mu_prior_xy + sqrt(var_prior)*randn()   ;
%             x(3*n  ,1)       = mu_prior_z  + sqrt(var_prior)*randn()   ;
         end
     end

     for n=1:Number_particles
         if b(n)==0
            for i = 2:sign_siz
                kxy            = sqrt(2*D*Delta_t(i))/rc           ;
                kz            = sqrt(2*D*Delta_t(i))/zc           ;
                
                x(n,i)         = x(n,i-1) + kxy*randn()           ;
                y(n,i)         = y(n,i-1) + kxy*randn()           ;
                z(n,i)         = z(n,i-1) + kz*randn()           ;
%                 x(3*n-1,i)         = x(3*n-1,i-1) + kxy*randn()                         ;
%                 x(3*n,i)           = x(3*n,i-1)   + kz *randn()                         ;
             end
        end
     end

end
