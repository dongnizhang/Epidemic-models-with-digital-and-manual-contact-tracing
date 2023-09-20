function counts_inf=epi_manual_ct(beta, gamma, delta, p, n, S0,I0,R0,D0,tf)


%%%%%%%%%%%%%%%%% Model Parameters:$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

time = 0; %time of event
S=S0; % number of susceptibles
I=I0;% number of infectives
R=R0;% number of recovered
D=D0;% number of diagnosed

% counts for those who have been infected
counts_inf=I0;
% counts for diagnosis due to themselves
counts_dia=0;
% counts for diagnosis due to being traced:
counts_diact=0;
% counts for the contact tracing events:
counts_ct=0;  

% construct a table for the number of S,I,R,D at each time of event t: 
% time array of events 
ta=cell(1,1);
ta{1,1}=time;
% the number of S at time
S_t=cell(1,1);
S_t{1,1}=S0;
% the number of I at time
I_t=cell(1,1);
I_t{1,1}=I0;
% the number of R at time
R_t=cell(1,1);
R_t{1,1}=R0;
% the number of D at time
D_t=cell(1,1);
D_t{1,1}=D0;

%label of Susceptibles
labelS=cell(S0,1);
for i=1:S0
    labelS{i,1}=i;
end
% label of Infectives 
labelI=cell(I0,1);
for i=1:I0
    labelI{i,1}=S0+i;
end
% label of Recovered
labelR=cell(R0,1);
for i=1:R0
    labelR{i,1}=S0+I0+i;
end
% label of Diagnosed
labelD=cell(D0,1);
for i=1:D0
    labelD{i,1}=S0+I0+R0+i;
end

% construct a table of contacts and sym-report contacts:

% the column for labeled individuals:
individual=(1:n)';
% the list of contacts of each individual: 
contacts=cell(n,1);
% the list of symmetric reporting contacts of each individual: 
symcon=cell(n,1);

T=table(individual,contacts,symcon);


% start the loop
while time < tf 
% First, we set the rates: 
    %fprintf('\n previous time %d\n',time);
% New infection:
    r_I=beta*I*(S/n); 
% Natural recovery:
    r_N=gamma*I;
% Diagnosis:    
    r_D=delta*I; 
% total rates:    
    r=r_I+r_N+r_D;
    
% time until next event happens:

    dt=exprnd(1/r);
    
    time=time+dt;
    
    %fprintf('\n current time %d\n',time);
    
    ta{size(ta,1)+1,1}=time;
    
if time >= tf
    break
end
        
%% Now, randomly choose a uniform(0,1) number to decide which event exactly occurs
    
    rd=rand;
    %fprintf('\n event random number %d\n',rd);
    %fprintf('\n New infection %d\n',r_I/r);
    %fprintf('\n Natural recovery %d\n',(r_I+r_N)/r);
%% --------------------------------First, if new infection occurs--------------------------------------
if rd <= (r_I/r) 
       
%% Randomly choose one of the remaining susceptibles to be infected:

       % random choose one position from the labeled susceptibles
       ind1=randperm(size(labelS,1),1);
       % get exactly which susceptible it is
       ri1=labelS{ind1,1};
       %fprintf('\n which S got infected %d\n',ri1);
       % delete this susceptible from the list of susceptibles
       labelS(ind1,:)=[];
       % record the number of S,I,R,D at this time of event t:  
       S_t{size(S_t,1)+1,1}=S-1;
       I_t{size(I_t,1)+1,1}=I+1;
       R_t{size(R_t,1)+1,1}=R;
       D_t{size(D_t,1)+1,1}=D;
       % the number of suscep. decreased by 1
       S=S-1;
       %fprintf('\n the number of S %d\n',S);
%% Randomly choose one of the infectives as the corr. infector 
      
       % random pick one of the infetives as the one who infected the ri1
        ind2=randperm(size(labelI,1),1);
       % get exactly which infective it is
        ri2=labelI{ind2,1};
        %fprintf('\n the ri1 got infected by this I %d\n',ri2);
       % i.e. the ri1 is infected by ri2. 
       
       % Add this new infected ri1 to the labelI 
        labelI{size(labelI,1)+1,1}=ri1;
        counts_inf=counts_inf+1;
        I=I+1;
       % Insert the contact between ri1 as infectee and ri2 as infector
        contacts{ri1,1}=[contacts{ri1,1},ri2];
        contacts{ri2,1}=[ contacts{ri2,1},ri1];
       
%% Decide whether ri1 reports ri2 and whether ri2 reports ri1 (Symmetric Reporting):
        
       % Random choose a uniform(0,1) number
           rdm = rand;
           %fprintf('\n reporting random number %d\n',rdm);
           %fprintf('\n the reporting probability %d\n',p);
       % If random smaller than p, then ri1(infectee) reports ri2(infector)
           if rdm<p
       % In the symcontacts table,
       % insert ri2(infector) under ri1(infectee)'s contact list
       % every infector is at the top of the corr. contact list 
            symcon{ri1,1}=[ri2,symcon{ri1,1}];
       %  insert ri1 under ri2's contact list
            symcon{ri2,1}=[symcon{ri2,1},ri1];
           else
               
           end
       T=table(individual,contacts,symcon);
           
           if I==0 %%%%%%%%%if there is no infectives%%%%%%%%%%%%
       
           break
           end
        
%% ------------------------------------if Natural Recovery occurs ----------------------------------
    elseif rd > (r_I/r) && rd <= ((r_I+r_N)/r)
%%  Randomly choose one of the remaining infectives to be naturally recovered:

% random choose one from the labeled infectives
       ind3=randperm(size(labelI,1),1);
% get exactly which infective it is 
       ri3=labelI{ind3,1}; 
       %fprintf('\n which I is naturally recovered %d\n',ri3);
% record the number of S,I,R,D at this time of event t:  
       S_t{size(S_t,1)+1,1}=S;
       I_t{size(I_t,1)+1,1}=I-1;
       R_t{size(R_t,1)+1,1}=R+1;
       D_t{size(D_t,1)+1,1}=D;
% So, the number of infectives decreased by 1. 
       I=I-1;
       %fprintf('\n the number of infectives %d\n',I);
% So, The number of recovered is increased by 1.
       R=R+1; 
       %fprintf('\n the number of recovered %d\n',R);
% Remove this recovered(no more infectious) from the lablelI: 
       labelI(ind3,:)=[];
% Add this new recovered to the labelR.
       labelR{size(labelR,1)+1,1}=ri3;
       
        if I==0 %%%%%%%%%if there is no infectives%%%%%%%%%%%%
       
        break
        end 
       
 
%% if Diagnosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    else
%%  Randomly choose one of the remaining infectives to be diagnosed:

% Random choose one from the labeled infectives
      ind4=randperm(size(labelI,1),1);
      
% get exactly which infective it is      
      ri4=labelI{ind4,1};
      
      %fprintf('\n which I is diagnosed %d\n',ri4);

       
% Remove this diagnosed(no more infectious) from the lablelI 
      labelI(ind4,:)=[];
      
% Add this new diagnosed to the labelD.
      labelD{size(labelD,1)+1,1}=ri4;
      
% the diagnosis is due to himself: 

      counts_dia=counts_dia+1;
      
%% ------------- symmetric Reporting Contacts Tracing ------------------

Diag=[ri4];% construct a list of individuals whose contacts will be traced.
counts=0 ; % this only counts Diagnosis due to Contact Tracing only in this time        
OldDiag=Diag;
% while the list is not empty:
while ~isempty(Diag) 
           
%fprintf('\n Diagosed list %d\n',Diag);   
    
              if  isempty(labelI) % if no infectives, stop the loop
                      break
              end
    
             NewDiag=[]; % ready to construct a new list of individuals whose contacts will be traced.
       
%------for each individual i from list Diag:
       
       for i = 1:length(Diag)
           
           if  isempty(labelI) % if no infectives, stop the loop
                      break
           end
         
           %fprintf('\nIteration i=%d\n',i);
           
           %fprintf('\n i-th individual will be traced %d\n',Diag(i));
           
           %fprintf('\n the symcontacts of i-th individual %d\n',symcon{Diag(i),1});
           
%------- if this i-th individual has symmetric contacts, then:   
           %fprintf('\n Does this i-th individual have symmetric contacts? %d\n',~isempty(symcon{Diag(i),1}));
           
           if ~isempty(symcon{Diag(i),1})
               
%------- for each j-th contact of the  i-th individual :

               for j=1:length(symcon{Diag(i),1})
                
                   %fprintf('\nIteration j=%d\n',j);
                   %fprintf('\nj-th contact of the  i-th individual is %d\n',symcon{Diag(i),1}(j));
                
                if  isempty(labelI) % if no infectives, stop the loop
                      break
                end
                
 %fprintf('\n Does this j-th contact have NOT been already diagnosed? %d\n',isempty(find([labelD{:}] == symcon{Diag(i),1}(j))));
               
%-------- if the j-th contact has not been already diagnosed, then continue it with next follow steps: 

               if isempty(find([labelD{:}] == symcon{Diag(i),1}(j)))
                  
%-------- then there can be two cases: either it is still infectious or it is already naturally recovered: 
             
%-------- if the j-th contact is still infectious, then:

               if  ~isempty(find([labelI{:}] == symcon{Diag(i),1}(j)))
                  
                  %fprintf('\n the j contact is still infectious? %d\n',~isempty(find([labelI{:}] == symcon{Diag(i),1}(j))));
                  
                  % find out at which position of the list labelI:
                  
                   i_j=find([labelI{:}] == symcon{Diag(i),1}(j));
                  
                  % Remove this j-th contact from infectives list labelI:
                  
                  labelI(i_j,:)=[];
                  
                  % Add this j-th contact to the diagnosed list:
                  
                  labelD{size(labelD,1)+1,1}=symcon{Diag(i),1}(j);
                  
                  % there is a Diagnosis event due to CT:
                  
                  counts=counts+1;
                  
                  counts_diact=counts_diact+1;
                  
                  
                  % Add this to be the new list NewDiag
                  
                  NewDiag=[NewDiag, symcon{Diag(i),1}(j)];
                  
                  % there is one contact tracing event occured. 
                  
                  %fprintf('\n counts_ct %d\n',counts_ct);
                  
                  counts_ct=counts_ct+1;
               
%-------- otherwise, this j-th contact is already naturally recovered, and we still tracing own contacts.
              
               else
 %fprintf('\n Does this j-th contact have been recovered? %d\n',~isempty(find([labelR{:}] == symcon{Diag(i),1}(j))));                 
                  %Since we still tracing its own contacts, add this to 
                  % the new list
                  
                  NewDiag=[NewDiag, symcon{Diag(i),1}(j)];
                  
                  % there is one CT event occured. 
                   
                  counts_ct=counts_ct+1;
                  
               
               

               end %-------- Row 272  

               else %-------- Row 266 Otherwise, the jth contact has been diagnosed, the tracing of this contac stops.
                   
               end %-------- Row 266
               
     
               end %-------- Row 254
               
           end %-------- Row 250
           
               
       end %-------- Row 236: the tracing of all from the list Diag stops
           
 %-------- Now, we compare the new list NewDiag and the old list OldDiag:
       
       % if some elements in the NewDiag are also in the OldDiag, then
       % there is no need to trace them again from beginning, so delete
       % them from the NewDiag:
       
       %fprintf('\n the old list %d\n',OldDiag);
       %fprintf('\n the new list %d\n',NewDiag);
       
       if ~isempty(intersect(OldDiag,NewDiag))
       
       % get the same elements from NewDiag and OldDiag:
       
       inter=intersect(OldDiag,NewDiag);
       
       % delete them from NewDiag
       
       for i=1:length(inter)
           indice=find(NewDiag==inter(i));
           NewDiag(indice)=[];
          
       end
       
       %fprintf('\n the cleared new list %d\n',NewDiag);
      
       else % otherwise, there is no intersection between them
       
       % we add these to the OldDiag, in order to compare for next time    
       OldDiag=[OldDiag,NewDiag];
       %fprintf('\n the extended old list %d\n',OldDiag);
       
       end %-------- Row 341
       
%-------- Now, we can restart the whole loop of tracing of this new list: 

       Diag = NewDiag;
      
       %fprintf('\n now we start tracing the updated list %d\n',Diag);
       
end %-------- Row 224: This is end for the whole tracing at this event time t.

%-------- In conclusion, Record the number of S,I,R,D at this time of event t:

       S_t{size(S_t,1)+1,1}=S;
       I_t{size(I_t,1)+1,1}=I-counts-1;
       R_t{size(R_t,1)+1,1}=R;
       D_t{size(D_t,1)+1,1}=D+counts+1;
     
      % the accumulated number of infectives:
      I=I-counts-1;
      
      % the accumulated number of diagnosed:
      D=D+counts+1;
      
  
  
      if I==0 %%%%%%%%%if there is no infectives, we exit the loop %%%%%%%%%%%%
       
        break
      end
      
      %toc
end


end

end
