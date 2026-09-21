% User Guidline:

% FFBS          ------> Forward Filtering Backward Sampling
% D             ------> Diffusion coefficent
% Lambda_P      ------> rate of photons that emitted by particle
% Lambda_bcgd   ------> Intensity or number of photons in the background


function [ D        , x                       , y                 , z                        , ...
           Lambda_P , acceptance_rate_lambdaP , Lambda_b          , acceptance_rate_lambda_b , ...
           b        , q                       , acceptance_rate_b , concentration          ] = ...
...
...
Gibbs_sampling_function_withNoise_Continous2(  ...
...
...
DeltaT               ,  Max_iter             ,  Prior_D_alpha    ,  Prior_D_beta         , ...
D                    ,  Lambda_P             ,  Lambda_b         ,  rc                   , ...
zc                   ,  Number_species       ,  alpha_q          ,  beta_q               , ...
b                    ,  q                    ,  M0_xy            ,  var0_xy              , ...
M0_z                 ,  var0_z               , ...
Prior_Lambda_P_alpha ,  Prior_Lambda_P_beta  ,  Lambda_P_alpha   ,  Prior_Lambda_b_alpha , ...
Prior_Lambda_b_beta  ,  Lambda_b_alpha       ,  x                ,  y                    , ...
z                    ,  concentration        ,  concen_radious       )


    rng('shuffle')
    sign_siz   = length(DeltaT); % size of the signal

    
    acceptance_rate_lambdaP  = 0  ;
    acceptance_rate_lambda_b = 0  ;
    acceptance_rate_b        = 0  ;
    

% Starting Loop for Gibbs sampling iteration for X , D of each one of the species
    for i=length(D)+1:length(D)+Max_iter

        for ki=1:Number_species 
%-------------- functions for learning Active locations X coordinates--------- 
           if b(ki,i-1)==1      %FFBS from active particle
              bb=b(:,i-1);
              bb(ki)=0;
              % calculate the effect of the other particles on the likelihood
              constant1 =  Lambda_b(i-1)+Lambda_P(i-1)*(bb'*exp(-2*((x.^2)+(y.^2)+(z.^2)) )); 
               
              % Updating  X
              x(ki,:) = FFBS_X( constant1 , DeltaT(1,:) , D(i-1) , rc , zc , Lambda_P(i-1)  , Lambda_b(i-1) , sign_siz , M0_xy , var0_xy , y(ki,:) , z(ki,:));
              % Updating  Y
              y(ki,:) = FFBS_Y( constant1 , DeltaT(1,:) , D(i-1) , rc , zc , Lambda_P(i-1)  , Lambda_b(i-1) , sign_siz , M0_xy , var0_xy , x(ki,:) , z(ki,:));
              % Updating  Z
              z(ki,:) = FFBS_Z( constant1 , DeltaT(1,:) , D(i-1) , rc , zc , Lambda_P(i-1)  , Lambda_b(i-1) , sign_siz , M0_z  , var0_z  , x(ki,:) , y(ki,:));
           end
        end
        
    
        % Calculate the concentration of particles          
        for l=1:length(concen_radious)
                 concentration(i-(floor((i-1)/1000)*1000),:,l) = ...
                                      b(:,i-1)'*(sqrt(((x).^2)+((y).^2)+((z).^2)) <=concen_radious(l)) ;
        end
        

               % inactive particle 
         [x , y, z]  =    Metropolis_inactive_Diff ( x  , y ,  z  , sign_siz , Number_species  , ...
                 DeltaT      , b(:,i-1)   , M0_xy     ,var0_xy  ,M0_z     ,var0_z  , rc  , zc  , D(i-1)) ;

     
     % Updating  D
            D(i)    = MH_D( DeltaT(1,:) , x , y , z , Prior_D_alpha, Prior_D_beta, sign_siz, Number_species,rc,zc);
       
     % Defining hyper-parameter qn for the load
            
            Lambda_P(i) = Lambda_P(i-1);
            b(:,i)      = b(:,i-1);
           for kk=1:10 
     % Updating the loads
           q(:,i)  = betarnd(alpha_q/Number_species+b(:,i),(beta_q*(Number_species-1)/Number_species)+1-b(:,i));
           
           [ b(:,i) , acceptance_rate_b ] = ...
               ...
            Load_BP( b(:,i)        , q(:,i)             ,  Lambda_P(i)  ,  rc        , ...
                     zc              , Lambda_b(i-1)      ,  x              ,  y         , ...
                     z               , acceptance_rate_b  ,  DeltaT(1,:)    ,  sign_siz  , ...
                     Number_species ); 
       
     % Updating the emission rate
          
           [ Lambda_P(i)  , acceptance_rate_lambdaP ] =  MH_Lambda_P (b(:,i), DeltaT(1,:) , Lambda_P(i),  x , y , z , rc  , zc , Lambda_b(i-1)  , ...
                               Prior_Lambda_P_alpha , Prior_Lambda_P_beta , Lambda_P_alpha, sign_siz , acceptance_rate_lambdaP);
           end
     % Updating the emission background
         
%            [ Lambda_b(i)  , acceptance_rate_lambda_b ] = MH_Lambda_b (b(:,i), DeltaT(1,:) , Lambda_b(i-1) , x , y , z , rc , zc , ...
%                Prior_Lambda_b_alpha , Prior_Lambda_b_beta , Lambda_b_alpha, Lambda_P(i), sign_siz , acceptance_rate_lambda_b); 

              Lambda_b(i)= Lambda_b(i-1);
              
        
         
            
    end

%     Xss=nan(Number_species,sign_siz);
%     Xss=sum(abs(Xs(:,:,2:end)),3)/(Max_iter-1);
%     
%     Yss=nan(Number_species,sign_siz);
%     Yss=sum(abs(Xs(:,:,2:end)),3)/(Max_iter-1);
%     
%     Zss=nan(Number_species,sign_siz);
%     Zss=sum(abs(Xs(:,:,2:end)),3)/(Max_iter-1);
    
    acceptance_rate_lambdaP = acceptance_rate_lambdaP/(2*(Max_iter-1));
    acceptance_rate_lambda_b = acceptance_rate_lambda_b/(2*(Max_iter-1));
    acceptance_rate_b = acceptance_rate_b/(2*(Max_iter-1));
    
end



%% FFBS for the hidden states, X
% 
function [ x_sampled ] = FFBS_X(constant1 , DeltaT , D , rc , zc , lambda_P , lambda_b , len , M0 , var0 , y, z)

M1=nan(1,length(DeltaT));
C1=nan(1,length(DeltaT));


% mean and variance of our normal approximation, equation 5, from Ioannis' write up
mu0= 8/9;
sig0= 1/9;

D= D/(rc^2);
% initialization of both mean and variance of normal

% lambda_b = 10000;

    lambdaXM0  =  constant1(1) + lambda_P*exp(-2*((M0).^2 + (y(1)).^2 + (z(1)).^2))      ;
    DlambdaXM0 = -4*M0*lambda_P*exp(-2*((M0).^2 + (y(1)).^2 + (z(1)).^2))   ;
    f        =  ((DeltaT(1)*lambdaXM0))^(1/3)                                                   ;
    fp       =  (((DeltaT(1))^(1/3))*DlambdaXM0*(lambdaXM0^(-2/3))/3) + eps                     ;
    a1      =  DlambdaXM0/(3*lambdaXM0)                                                         ;
    b1       =  (mu0-f+M0*fp)/fp                                                                ;
% variance of normal ( according to write up c(n) )
    C1(1)    =  (1./((1./(var0))+((fp^2/sig0))));
% mean of normal ( according to write up m(n) )
    M1(1)    =  ( a1 + ((fp^2).*b1./sig0)  +  (M0./(var0)) )*C1(1);

for i=2:len
%--------------------- at the point m(n-1)
    lambdaX  =  constant1(i-1) + lambda_P*exp(-2*((M1(i-1)).^2 + (y(i)).^2 + (z(i)).^2))     ;
    DlambdaX = -4*M1(i-1)*lambda_P*exp(-2*((M1(i-1)).^2 + (y(i)).^2 + (z(i)).^2))   ;
    f        =  ((DeltaT(i)*lambdaX))^(1/3)                                                          ;
    fp       =  (((DeltaT(i))^(1/3))*DlambdaX*(lambdaX^(-2/3))/3) + eps                              ;
    CC1      =  DlambdaX/(3*lambdaX)                                                                 ;
    h1       =  (mu0-f+M1(i-1)*fp)/fp                                                                ;
% variance of normal ( according to write up c(n) )
    C1(i)    =  (1./((1./(C1(i-1)+(2*D.*DeltaT(i))))+((fp^2/sig0))));
% mean of normal ( according to write up m(n) )
    M1(i)    =  ( CC1 + ((fp^2).*h1./sig0)  +  (M1(i-1)./(C1(i-1)+(2*D.*DeltaT(i)))) )*C1(i);
end

%--------------------- Sampling from last filter An  directly

% BECAREFUL that the filter is symmetric relative to the point of origin and
% because of that we need to sample from sum of two normal in the backward step
     if  rand()>0.5
         x_sampled(len) =  M1(len)+sqrt(C1(len))*randn();                                    
     else
         x_sampled(len) = -M1(len)+sqrt(C1(len))*randn();        
     end
%----------------- Sampling from the rest of filters, An-1, An-2, ..., till A1
     for ii= len-1:-1:1  
         % normalization factors of each normal
         if  rand()< 1/( 1+exp( -2*x_sampled(ii+1)*M1(ii)/((C1(ii)+2*D.*DeltaT(ii+1))) ) )
             x_sampled(ii)= ((M1(ii)*2*D.*DeltaT(ii+1) + C1(ii).*x_sampled(ii+1))./(C1(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./C1(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         else           
             x_sampled(ii)= ((-M1(ii)*2*D.*DeltaT(ii+1) + C1(1,ii).*x_sampled(ii+1))./(C1(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./C1(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         end
     end            
end

%% FFBS for the hidden states, Y

function [ y_sampled ] = FFBS_Y(constant1 , DeltaT , D , rc , zc , lambda_P , lambda_b , len , M0 , var0 , x , z)

MY=nan(1,length(DeltaT));
CY=nan(1,length(DeltaT));
% mean and variance of our normal approximation, equation 5, from Ioannis' write up
mu0= 8/9;
sig0= 1/9;

D=(D/(rc^2));
% initialization of both mean and variance of normal

% lambda_b = 10000;

    lambdaYM0  =  constant1(1) + lambda_P*exp(-2*((M0).^2 + (x(1)).^2 + (z(1)).^2))      ;
    DlambdaYM0 = -4*M0*lambda_P*exp(-2*((M0).^2 + (x(1)).^2 + (z(1)).^2))   ;
    f        =  ((DeltaT(1)*lambdaYM0))^(1/3)                                                   ;
    fp       =  (((DeltaT(1))^(1/3))*DlambdaYM0*(lambdaYM0^(-2/3))/3) + eps                     ;
    a1      =  DlambdaYM0/(3*lambdaYM0)                                                         ;
    b1       =  (mu0-f+M0*fp)/fp                                                                ;
% variance of normal ( according to write up c(n) )
    CY(1)    =  (1./((1./(var0))+((fp^2/sig0))));
% mean of normal ( according to write up m(n) )
    MY(1)    =  ( a1 + ((fp^2).*b1./sig0)  +  (M0./(var0)) )*CY(1);

for i=2:len
%--------------------- at the point m(n-1)
lambdaY = constant1(i) + lambda_P*exp(-2*((MY(i-1))^2 + (x(i)).^2 + (z(i)).^2))        ;
DlambdaY = -4*MY(i-1)*lambda_P*exp(-2*((MY(i-1))^2 + (x(i)).^2 + (z(i)).^2))  ;   
f= ((DeltaT(i)*lambdaY))^(1/3)                                                              ;     
fp= (((DeltaT(i))^(1/3))*DlambdaY*(lambdaY^(-2/3))/3) + eps                                 ; 
CC1Y = DlambdaY/(3*lambdaY)                                                                 ;
h1 =(mu0-f+MY(i-1)*fp)/fp                                                                   ;
% variance of normal ( according to write up c(n) )
CY(i) = (1./((1./(CY(i-1)+(2*D.*DeltaT(i))))+((fp^2/sig0))))                              ;
% mean of normal ( according to write up m(n) )
MY(i) = ( CC1Y + ((fp^2).*h1./sig0)  +  (MY(i-1)./(CY(i-1)+(2*D.*DeltaT(i)))) )*CY(i);
end

  
%--------------------- Sampling from last filter An  directly

% BECAREFUL that the filter is symmetric relative to the point of origin 
     if  rand()>0.5
         y_sampled(len) =  MY(len)+sqrt(CY(len))*randn();                                    
     else
         y_sampled(len) = -MY(len)+sqrt(CY(len))*randn();        
     end
%----------------- Sampling from the rest of filters, An-1, An-2, ..., till A1
     for ii= len-1:-1:1  
         % normalization factors of each normal
         if  rand()< 1/( 1+exp( -2*y_sampled(ii+1)*MY(ii)/((CY(ii)+2*D.*DeltaT(ii+1))) ) )
             y_sampled(ii)= ((MY(ii)*2*D.*DeltaT(ii+1) + CY(ii).*y_sampled(ii+1))./(CY(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./CY(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         else           
             y_sampled(ii)= ((-MY(ii)*2*D.*DeltaT(ii+1) + CY(1,ii).*y_sampled(ii+1))./(CY(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./CY(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         end
     end 
end


%% FFBS for the hidden states, Z
% 
function [ z_sampled ] = FFBS_Z(constant1 , DeltaT , D , rc , zc , lambda_P , lambda_b  , len , M0 , var0 , x , y )

MZ=nan(1,length(DeltaT));
CZ=nan(1,length(DeltaT));
% mean and variance of our normal approximation, equation 5, from Ioannis' write up
mu0= 8/9;
sig0= 1/9;


D=(D/(zc^2));

% initialization of both mean and variance of normal
% lambda_b = 10000;

    lambdaZM0  =  constant1(1) + lambda_P*exp(-2*((M0).^2 + (x(1)).^2 + (y(1)).^2))      ;
    DlambdaZM0 = -4*M0*lambda_P*exp(-2*((M0).^2 + (x(1)).^2 + (y(1)).^2))   ;
    f        =  ((DeltaT(1)*lambdaZM0))^(1/3)                                                   ;
    fp       =  (((DeltaT(1))^(1/3))*DlambdaZM0*(lambdaZM0^(-2/3))/3) + eps                     ;
    a1      =  DlambdaZM0/(3*lambdaZM0)                                                         ;
    b1       =  (mu0-f+M0*fp)/fp                                                                ;
% variance of normal ( according to write up c(n) )
    CZ(1)    =  (1./((1./(var0))+((fp^2/sig0))));
% mean of normal ( according to write up m(n) )
    MZ(1)    =  ( a1 + ((fp^2).*b1./sig0)  +  (M0./(var0)) )*CZ(1);

for i=2:len
%--------------------- at the point m(n-1)
lambdaZ = constant1(i)  + lambda_P*exp(-2*((MZ(i-1))^2 + (x(i)).^2 + (y(i)).^2))       ; 
DlambdaZ = -4*MZ(i-1)*lambda_P*exp(-2*((MZ(i-1))^2 + (x(i)).^2 + (y(i)).^2))  ;    
f= ((DeltaT(i)*lambdaZ))^(1/3)                                                              ;
fp= (((DeltaT(i))^(1/3))*DlambdaZ*(lambdaZ^(-2/3))/3) + eps                                 ; 
CC1Z = DlambdaZ/(3*lambdaZ)                                                                      ;
h1 =(mu0-f+MZ(i-1)*fp)/fp                                                                       ;
% variance of normal ( according to write up c(n) )
CZ(i) = (1./((1./(CZ(i-1)+(2*D.*DeltaT(i))))+((fp^2/sig0))))                              ;
% mean of normal ( according to write up m(n) )
MZ(i) = ( CC1Z + ((fp^2).*h1./sig0)  +  (MZ(i-1)./(CZ(i-1)+(2*D.*DeltaT(i)))) )*CZ(i)      ;
end

%--------------------- Sampling from last filter An  directly
% BECAREFUL that the filter is symmetric relative to the point of origin
     if  rand()>0.5
         z_sampled(len) =  MZ(len)+sqrt(CZ(len))*randn();                                    
     else
         z_sampled(len) = -MZ(len)+sqrt(CZ(len))*randn();        
     end
%----------------- Sampling from the rest of filters, An-1, An-2, ..., till A1
     for ii= len-1:-1:1  
         % normalization factors of each normal
         if  rand()< 1/( 1+exp( -2*z_sampled(ii+1)*MZ(ii)/((CZ(ii)+2*D.*DeltaT(ii+1))) ) )
             z_sampled(ii)= ((MZ(ii)*2*D.*DeltaT(ii+1) + CZ(ii).*z_sampled(ii+1))./(CZ(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./CZ(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         else           
             z_sampled(ii)= ((-MZ(ii)*2*D.*DeltaT(ii+1) + CZ(1,ii).*z_sampled(ii+1))./(CZ(ii)+2*D.*DeltaT(ii+1)))+...
                             sqrt(1./((1./CZ(ii))+(1./(2*D.*DeltaT(ii+1)))))*randn();
         end
     end 
end


%% M for the D for X

function [ D_sampled ] = MH_D( DeltaT, x , y , z , Prior_D_alpha, Prior_D_beta, sign_siz, Number_species, rc,  zc)


    alpha1 = Prior_D_alpha + 1.5*Number_species*(sign_siz-2);
    beta1= 1./(Prior_D_beta+(sum( sum(  (rc^2)*((x(:,2:end-1)-x(:,3:end)).^2)+...
                                          (rc^2)*((y(:,2:end-1)-y(:,3:end)).^2)+...
                                          (zc^2)*((z(:,2:end-1)-z(:,3:end)).^2) ,1)./(DeltaT(3:end)) )/4));

                                      
    D_sampled = 1./gamrnd(alpha1,beta1);
end


%% MH for the rate lambda_P
function [ lambda_P_sampled, acceptance_rate_P] = MH_Lambda_P( b, DeltaT, lambda_P_old , x , y , z , rc , zc , lambda_b , Prior_lambda_P_alpha, Prior_lambda_P_beta, alpha,  sign_siz , acceptance_rate_P)
          
% keyboard
         lambda_P_new = gamrnd(alpha,lambda_P_old/alpha);
        
        for i=1:sign_siz-1
            
            pp(:,i) =sum(b.*exp(-2*((x(:,i)).^2 + (y(:,i)).^2 + (z(:,i)).^2)));
            pp1(:,i)=sum(b.*exp(-2*((x(:,i)).^2 + (y(:,i)).^2 + (z(:,i)).^2))).*DeltaT(:,i);
            pp2(:,i)=log((lambda_b+lambda_P_new*pp(:,i))./(lambda_b+lambda_P_old*pp(:,i)));


        end  
        
        logr = (sum(pp*(lambda_P_old-lambda_P_new)*DeltaT(:,i)))+ (((Prior_lambda_P_beta)^-1)*(lambda_P_old-lambda_P_new))+...
                sum(pp2)+(Prior_lambda_P_alpha-1)*log(lambda_P_new/lambda_P_old)...
                ...
               + ((2*alpha-1)*log(lambda_P_old/lambda_P_new)) + (alpha*((lambda_P_new/lambda_P_old)-(lambda_P_old/lambda_P_new))); 
           
 
         if  logr>(-exprnd(1))
             lambda_P_sampled = lambda_P_new;
             acceptance_rate_P=acceptance_rate_P+1;
         else
             lambda_P_sampled = lambda_P_old;
         end
end

%% MH for the rate lambda_b

function [ lambda_b_sampled, acceptance_rate_b] = MH_Lambda_b(b, DeltaT, lambda_b_old , x , y , z , rc , zc , Prior_lambda_b_alpha, Prior_lambda_b_beta, alpha, lambda_P, sign_siz , acceptance_rate_b)
          

         lambda_b_new = gamrnd(alpha,lambda_b_old/alpha);
        
        for i=1:sign_siz-1
            
            bb(:,i)=sum(b.*exp(-2*((x(:,i)).^2 + (y(:,i)).^2 + (z(:,i)).^2)));
%             bb1(:,i) = DeltaT(:,i);
            bb2(:,i)=log((lambda_b_new+lambda_P*bb(:,i))./(lambda_b_old+lambda_P*bb(:,i)));

        end  
        
        logr = ((sum(DeltaT))*(lambda_b_old-lambda_b_new))+ ((Prior_lambda_b_beta)^-1)*(lambda_b_old-lambda_b_new) + ...
            sum(bb2)+(Prior_lambda_b_alpha-1)*log(lambda_b_new/lambda_b_old)...
            + ((2*alpha-1)*log(lambda_b_old/lambda_b_new)) + (alpha*((lambda_b_new/lambda_b_old)-(lambda_b_old/lambda_b_new))); 
            
   
         if  logr>(-exprnd(1))
             lambda_b_sampled = lambda_b_new;
             acceptance_rate_b=acceptance_rate_b+1;
         else
             lambda_b_sampled = lambda_b_old;
         end
end
