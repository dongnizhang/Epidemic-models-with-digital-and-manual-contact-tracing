
rng(123);

% size of population: 

n=1000; %n=5000

% numbr of simulations:

n_sim=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epidemic without any tracing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the number of individuals have been infected: 

Z1 = cell(1,1);

% the major outbreak parts:

Z1_major=cell(1,1);

% the minor outbreak parts:

Z1_minor=cell(1,1);


% infection rate: beta
% rate of natural recovery: gamma
% rate of diagnosis: delta
% probability of manual contact tracing: p
% size of population: n
% number of initial susceptibles: S0
% number of initial infectives: I0 
% Note! I0 should be positive!
% number of initial recovered: R0
% number of initial diagnosed: D0
% time interval: [0, tf]

for i=1:n_sim
    
    Z1{i,1}=epi_without_ct(0.8,1/7,1/7,n,(n-1),1,0,0,100);
    
    if Z1{i,1} <= (0.1*n) % minor outbreak
       
        Z1_minor{size(Z1_minor,1)+1,1}=Z1{i,1};
        
        
    else % major outbreak
        
        Z1_major{size(Z1_major,1)+1,1}=Z1{i,1};
        
    end
     
end





%%I_sym_minor{1,1}=[];
%% sample mean fraction for major outbreaks

m1_sym=mean(cell2mat(Z1_major)/n);
disp('the mean fraction without ct')
disp(m1_sym)

%% standard deviation for major outbreaks

s1_sym=std(cell2mat(Z1_major)/n);
disp('the standard deviation without ct')
disp(s1_sym)

%% the fraction of minor/major outbreak:

disp('the fraction of minor outbreaks without ct')
disp((size(Z1_minor,1)-1)/10000);

disp('the fraction of major outbreaks without ct')
disp((size(Z1_major,1)-1)/10000);

%% display the histogram of all simulations:

h=histogram(cell2mat(Z1),'BinWidth',10);
ylim([0,10000]);
xlim([0, n]);
%saveas(h1, 'Histo_sym.jpg', 'jpg');
%'FaceAlpha',1

h_major=histogram(cell2mat(Z1_major),'BinWidth',0.1*std(cell2mat(Z1_major)));
%ylim([0,200]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epidemic with manual tracing only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the number of individuals have been infected: 

Z2 = cell(1,1);

% the major outbreak parts:

Z2_major=cell(1,1);

% the minor outbreak parts:

Z2_minor=cell(1,1);

% infection rate: beta
% rate of natural recovery: gamma
% rate of diagnosis: delta
% probability of manual contact tracing: p
% size of population: n
% number of initial susceptibles: S0
% number of initial infectives: I0 
% Note! I0 should be positive!
% number of initial recovered: R0
% number of initial diagnosed: D0
% time interval: [0, tf]

for i=1:n_sim
    
    Z2{i,1}=epi_manual_ct(0.8,1/7,1/7,2/3,n,(n-1),1,0,0,100);
   
    %disp(i)
    
    if Z2{i,1} <= (0.1*n) % minor outbreak
       
        Z2_minor{size(Z2_minor,1)+1,1}=Z2{i,1};
        
        
    else % major outbreak
        
        Z2_major{size(Z2_major,1)+1,1}=Z2{i,1};
        
    end
     
end





%%I_sym_minor{1,1}=[];
%% sample mean fraction for major outbreaks

m2_sym=mean(cell2mat(Z2_major)/n);
disp('the mean fraction for manual ct')
disp(m2_sym)

%% standard deviation for major outbreaks

s2_sym=std(cell2mat(Z2_major)/n);
disp('the standard deviation for manual ct')
disp(s2_sym)

%% the fraction of minor/major outbreak:

disp('the fraction of minor outbreaks for manual ct')
disp((size(Z2_minor,1)-1)/10000);

disp('the fraction of major outbreaks for manual ct')
disp((size(Z2_major,1)-1)/10000);

%% display the histogram of all simulations:

h_manual=histogram(cell2mat(Z2),'BinWidth',10);
ylim([0,10000]);
xlim([0, n]);
%'FaceAlpha',1

%% display the histogram of all simulations:
h_manual_major=histogram(cell2mat(Z2_major),'BinWidth',0.1*std(cell2mat(Z2_major)));
%ylim([0,200]); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epidemic with digital tracing only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the number of individuals have been infected: 

Z3 = cell(1,1);

% the major outbreak parts:

Z3_major=cell(1,1);

% the minor outbreak parts:

Z3_minor=cell(1,1);

% infection rate: beta
% rate of natural recovery: gamma
% rate of diagnosis: delta
% probability of digital contact tracing: p_A
% size of population: n
% number of initial susceptibles: S0
% number of initial infectives: I0 
% Note! I0 should be positive!
% number of initial recovered: R0
% number of initial diagnosed: D0
% time interval: [0, tf]

for i=1:n_sim
    %disp(i)
    Z3{i,1}=epi_digital_ct(0.8,1/7,1/7,2/3,n,(n-1),1,0,0);
    
    if Z3{i,1} <= (0.1*n) % minor outbreak
       
        Z3_minor{size(Z3_minor,1)+1,1}=Z3{i,1};
        
        
    else % major outbreak
        
        Z3_major{size(Z3_major,1)+1,1}=Z3{i,1};
        
    end
     
end


%% sample mean fraction for major outbreaks

m3_sym=mean(cell2mat(Z3_major)/n);
disp('the mean fraction for digital ct')
disp(m3_sym)

%% standard deviation for major outbreaks

s3_sym=std(cell2mat(Z3_major)/n);
disp('the standard deviation for digital ct')
disp(s3_sym)

%% the fraction of minor/major outbreak:

disp('the fraction of minor outbreaks for digital ct')
disp((size(Z3_minor,1)-1)/10000);

disp('the fraction of major outbreaks for digital ct')
disp((size(Z3_major,1)-1)/10000);

%% display the histogram of all simulations:

h_digital=histogram(cell2mat(Z3),'BinWidth',10);
ylim([0,10000]);
xlim([0, n]);
%saveas(h1, 'Histo_sym.jpg', 'jpg');
%'FaceAlpha',1

%% display the histogram of all simulations:
h_digital_major=histogram(cell2mat(Z3_major),'BinWidth',0.1*std(cell2mat(Z3_major)));
%ylim([0,200]); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epidemic with manual ad digital tracing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the number of individuals have been infected: 

Z = cell(1,1);

% the major outbreak parts:

Z_major=cell(1,1);

% the minor outbreak parts:

Z_minor=cell(1,1);

% infection rate: beta
% rate of natural recovery: gamma
% rate of diagnosis: delta
% probability of digital contact tracing: p_A
% probability of manual contact tracing: p
% size of population: n
% number of initial susceptibles: S0
% number of initial infectives: I0 
% Note! I0 should be positive!
% number of initial recovered: R0
% number of initial diagnosed: D0
% time interval: [0, tf]

for i=1:n_sim
    
    Z{i,1}=epi_comb_ct(0.8,1/7,1/7,2/3,2/3,n,(n-1),1,0,0);
    
    if Z{i,1} <= (0.1*n) % minor outbreak
       
        Z_minor{size(Z_minor,1)+1,1}=Z{i,1};
        
        
    else % major outbreak
        
        Z_major{size(Z_major,1)+1,1}=Z{i,1};
        
    end
     
end


%% sample mean fraction for major outbreaks

m_sym=mean(cell2mat(Z_major)/n);
disp('the mean fraction  for combined ct')
disp(m_sym)

%% standard deviation for major outbreaks

s_sym=std(cell2mat(Z_major)/n);
disp('the standard deviation  for combined ct')
disp(s_sym)

%% the fraction of minor/major outbreak:

disp('the fraction of minor outbreaks for combined ct')
disp((size(Z_minor,1)-1)/10000);

disp('the fraction of major outbreaks for combined ct')
disp((size(Z_major,1)-1)/10000);

%% display the histogram of all simulations:

h_comb=histogram(cell2mat(Z),'BinWidth',10);
ylim([0,10000]);
xlim([0, n]);
%saveas(h1, 'Histo_sym.jpg', 'jpg');
%'FaceAlpha',1

%% display the histogram of all simulations:
h_comb_major=histogram(cell2mat(Z_major),'BinWidth',0.1*std(cell2mat(Z_major)));
%ylim([0,200]); 












% % theoretical mean: determinstic
% mu=mean(cell2mat(Z_major));
% 
% % standard deviation
% s=0;
% I_major=cell2mat(Z_major);
% n_major=size(I_major);
% n_major=n_major(1);
% for i=1:n_major
%     s=s+(I_major(i)-mu).^2;
% end
% sd=sqrt(s/(n_major-1));
% x_low=max(mu-3*sd,0);
% x_upper=mu+3*sd;
% xlim([x_low, x_upper]);
% hold on
% dh = h2.BinWidth;
% lo = min(I_major);
% hi = max(I_major);
% dxx = (hi-lo)/100;
% xx = linspace(lo,hi,101);
% pdf = 1/(2*pi*sd)*exp(-(xx-mu).^2/(2*sd^2));
% scalefactor = sum(h2.Values * dh)/(trapz(pdf)* dxx);
% pdf = scalefactor*pdf;
% plot(xx,pdf,'LineWidth',2)        
% hold off
% 
% histfit(cell2mat(Z_major),120)
% ylim([0,200]);
% x_low=max(mu-3*sd,0);
% x_upper=mu+3*sd;
% xlim([x_low, x_upper]);
% %pd_normal = fitdist(cell2mat(I_sym_major),'Normal')