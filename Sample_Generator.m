function[ Time , emission_rate , DeltaT , X0, XX , Y0 , YY , Z0 , ZZ] =  ...
...
...
Sample_Generator(minn, maxx, D, Lambda_P_real , Lambda_b_real ,  Rx , Ry , Rz , rc , zc , ...
Number_particles , Length_signal)
 
 

Time=zeros(Length_signal+1,1);
 
emission_rate=zeros(Length_signal+1,Number_particles);
 
x1 = zeros ( 2 , Number_particles );
y1 = zeros ( 2 , Number_particles );
z1 = zeros ( 2 , Number_particles );
 
X0 = nan (Number_particles,1);
Y0 = nan (Number_particles,1);
Z0 = nan (Number_particles,1);
 
XX = nan(Length_signal+1,Number_particles);
YY = nan(Length_signal+1,Number_particles);
ZZ = nan(Length_signal+1,Number_particles);
 
% Distribution of the particle randomly throughout of the space
 
for jk1=1:Number_particles
    
    axint1=4*Rx;
    ayint1=4*Ry;
    azint1=4*Rz;
    while (((axint1^2)+(ayint1^2)+(azint1^2)>(Rx^2 + Ry^2 + Rz^2)))
        axint1=(Rx-0.00000001)*(1-2*rand());
        ayint1=(Ry-0.00000001)*(1-2*rand());
        azint1=(Rz-0.00000001)*(1-2*rand());
    end
    x1(1,jk1) = axint1;
    y1(1,jk1) = ayint1;
    z1(1,jk1) = azint1;
    
    X0(jk1)=x1(1,jk1);
    XX(1,jk1)=X0(jk1);
    Y0(jk1)=y1(1,jk1);
    YY(1,jk1)=Y0(jk1);
    Z0(jk1)=z1(1,jk1);
    ZZ(1,jk1)=Z0(jk1);
end
 
 
 
 for i=2:Length_signal+1  
 
  
    emittion_rate_total(i-1)=sum(exp(-2*((x1(1,:)/rc).^2 + (y1(1,:)/rc).^2 + (z1(1,:)/zc).^2)));
    
    deltatt=exprnd(1/(Lambda_b_real+Lambda_P_real*emittion_rate_total(i-1)));
    
    for j1=1:Number_particles
        
        emission_rate(i-1,j1)=exp(-2*((x1(1,j1)/rc).^2 + (y1(1,j1)/rc).^2 + (z1(1,j1)/zc).^2));
        
        Time(i)=Time(i-1)+deltatt;
        % Brownian motion for the particle in each time step
        k1 =sqrt(2 * D * deltatt); 
        rkk1 = normrnd(x1(1,j1),k1);
        x1(2,j1)=rkk1(:,1);
        XX(i,j1)=x1(2,j1);
        
        rkk1 = normrnd(y1(1,j1),k1);
        y1(2,j1)=rkk1(:,1);
        YY(i,j1)=y1(2,j1);
        
        rkk1 = normrnd(z1(1,j1),k1);
        z1(2,j1)=rkk1(:,1);
        ZZ(i,j1)=z1(2,j1);
        
        % periodic boundaries
        if  x1(2,j1)>= Rx
            x1(2,j1)=x1(2,j1)-2*Rx;
            
        elseif x1(2,j1)<=-Rx
            x1(2,j1)=x1(2,j1)+2*Rx;
         
        end 
        if  y1(2,j1)>= Ry
            y1(2,j1)=y1(2,j1)-2*Ry;
            
        elseif y1(2,j1)<=-Ry
            y1(2,j1)=y1(2,j1)+2*Ry;
         
        end 
        
        if  z1(2,j1)>= Rz
            z1(2,j1)=z1(2,j1)-2*Rz;
            
        elseif z1(2,j1)<=-Rz
            z1(2,j1)=z1(2,j1)+2*Rz;
         
        end 
        
    end
 
    x1(1,:)=x1(2,:);
    y1(1,:)=y1(2,:);
    z1(1,:)=z1(2,:);
 end
 
SSignal_time=Time(min(find(Time>=minn))+find(Time(find(Time>=minn))<=maxx)-1);
DeltaT=(diff(SSignal_time))';
 
end

%%





























% function[ Time , emittion_rate , Delta_t , X0 , XX , Y0 , YY , Z0 , ZZ] =  ...
% ...
% ...
% Sample_Generator(D, Lambda_P_real , Lambda_b_real , R , rc , zc , Number_particles , Length_signal , time_end)
% 
% 
% Time=zeros(Length_signal+1,Number_particles);
%  
% emittion_rate=zeros(Length_signal+1,Number_particles);
% 
% x1 = zeros ( 2 , Number_particles );
% y1 = zeros ( 2 , Number_particles );
% z1 = zeros ( 2 , Number_particles );
% 
% X0 = nan (Number_particles,1);
% Y0 = nan (Number_particles,1);
% Z0 = nan (Number_particles,1);
% 
% XX = nan(Length_signal+1,Number_particles);
% YY = nan(Length_signal+1,Number_particles);
% ZZ = nan(Length_signal+1,Number_particles);
% 
% % Distribution of the particle randomly throughout of the space
% 
% for jk1=1:Number_particles
%     
%     axint1=4*R;
%     ayint1=4*R;
%     azint1=4*R;
%     while (((axint1^2)+(ayint1^2)+(azint1^2)>R^2))
%         axint1=(R-0.00000001)*(1-2*rand());
%         ayint1=(R-0.00000001)*(1-2*rand());
%         azint1=(R-0.00000001)*(1-2*rand());
%     end
%     x1(1,jk1) = axint1;
%     y1(1,jk1) = ayint1;
%     z1(1,jk1) = azint1;
%     
%     X0(jk1)=x1(1,jk1);
%     XX(1,jk1)=X0(jk1);
%     Y0(jk1)=y1(1,jk1);
%     YY(1,jk1)=Y0(jk1);
%     Z0(jk1)=z1(1,jk1);
%     ZZ(1,jk1)=Z0(jk1);
% end
%  
% 
%  
%  for i=2:Length_signal+1  
% 
%   
%     
%     for j1=1:Number_particles
%         
%         emittion_rate(i-1,j1)=exp(-2*((x1(1,j1)/rc).^2 + (y1(1,j1)/rc).^2 + (z1(1,j1)/zc).^2));
%         
%         deltatt=exprnd(1/(Lambda_P_real*emittion_rate(i-1,j1)));
%         Time(i,j1)=Time(i-1,j1)+deltatt;
%         % Brownian motion for the particle in each time step
%         k1 =sqrt(2 * D * deltatt); 
%         rkk1 = normrnd(x1(1,j1),k1);
%         x1(2,j1)=rkk1(:,1);
%         XX(i,j1)=x1(2,j1);
%         
%         rkk1 = normrnd(y1(1,j1),k1);
%         y1(2,j1)=rkk1(:,1);
%         YY(i,j1)=y1(2,j1);
%         
%         rkk1 = normrnd(z1(1,j1),k1);
%         z1(2,j1)=rkk1(:,1);
%         ZZ(i,j1)=z1(2,j1);
%         
%         % periodic boundaries
%         if  x1(2,j1)>= R
%             x1(2,j1)=x1(2,j1)-2*R;
%             
%         elseif x1(2,j1)<=-R
%             x1(2,j1)=x1(2,j1)+2*R;
%          
%         end 
%         if  y1(2,j1)>= R
%             y1(2,j1)=y1(2,j1)-2*R;
%             
%         elseif y1(2,j1)<=-R
%             y1(2,j1)=y1(2,j1)+2*R;
%          
%         end 
%         
%         if  z1(2,j1)>= R
%             z1(2,j1)=z1(2,j1)-2*R;
%             
%         elseif z1(2,j1)<=-R
%             z1(2,j1)=z1(2,j1)+2*R;
%          
%         end 
%         
%     end
%  
%     x1(1,:)=x1(2,:);
%     y1(1,:)=y1(2,:);
%     z1(1,:)=z1(2,:);
%  end
%  
%     for k=1:Number_particles
%         emittion_rate(end,k)=exp(-2*((x1(1,j1)/rc).^2 + (y1(1,j1)/rc).^2 + (z1(1,j1)/zc).^2));
%     end
%  
% % keyboard
%  
% Time_total=[];
% for i=1:Number_particles
%     Time_total= [Time_total , Time(:,i)'];
% end
% 
% Time_total=unique(sort(Time_total));
% 
% 
% j=1;
% time_back(1)=0;
% while time_back(j)<Time_total(end)
%       j=j+1;
%       deltatt=exprnd(1/Lambda_b_real);
%       time_back(1,j)=time_back(1,j-1)+deltatt;
% end
% 
% Time_total= [Time_total , time_back];
% 
% Time_total=unique(sort(Time_total));
% 
% 
% 
% jj=0;
% mn=0;
% for i=1:length(Time_total)
%     if Time_total(i)>time_end
%        if mn==0 
%           jj=i;
%           mn=1;
%        end
%     end
% end
%     
% Delta_t=diff(Time_total(1:jj));
% 
% end