function[ b , acceptance_rate_b ]    =...
...
Load_BP( ...
...
 b_old      ,   q   ,  Lambda_P    ,         rc           ,      zc     ,  Lambda_b  , ...
 x          ,   y   ,     z        ,  acceptance_rate_b   ,   Delta_t   ,  sign_siz  ,  Number_particles     ) 
 
% %-------- INPUT --------%
     % b_old                         :  Old load b
     % x                             :  Positions of the particles
     % acceptance_rate_b             :  Previous acceptance rate 
     % sign_siz                      :  Length of the observation signal
     % Delta_t                       :  Data acquesition time
     % Number_particles              :  Number of defined particles
     
%---------- OUTPUT -------%
     % b                             :  Updated load b
     % acceptance_rate_b             :  Updated acceptance rate
   
 
   
   % Set the vectors
   constant_new=zeros(1,sign_siz);
   constant_old=zeros(1,sign_siz);
   constant_old_new=zeros(1,sign_siz);
         


   % Propose the new loads
   b_new = rand(Number_particles,1)<q;
 
   % Calculate the constants  
   for n=1:Number_particles
       constant_old     = constant_old     + b_old(n)*exp(-2*((x(n,:)/rc).^2 + (y(n,:)/rc).^2 + (z(n,:)/zc).^2));
       constant_new     = constant_new     + b_new(n)*exp(-2*((x(n,:)/rc).^2 + (y(n,:)/rc).^2 + (z(n,:)/zc).^2));
       constant_old_new = constant_old_new + (b_old(n)-b_new(n))*exp(-2*((x(n,:)/rc).^2 + (y(n,:)/rc).^2 + (z(n,:)/zc).^2));
   end
 


   % Log of the Postrior ratio    
   logr= sum(log( (Lambda_b+Lambda_P*constant_new)./(Lambda_b+Lambda_P*constant_old)   )) + (Lambda_P*constant_old_new)*Delta_t';

   % Check the acceptance rate of load        
   if  logr>(-exprnd(1))
       b = b_new ;  
       acceptance_rate_b = acceptance_rate_b + 1 ;
   else
       b = b_old ;   
   end
   
end
