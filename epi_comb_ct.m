function counts_inf=epi_comb_ct(beta,gamma,delta,p_A,p,n,S0,I0,R0,D0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Model parameter:%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% n = size of population
%%% beta = transmission rate
%%% gamma = rate of natural recovery
%%% delta = rate of diagnosis(testing)
%%% p_A = probability of using tracing app
%%% p = probability of manual tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of initial susceptibles: S0
% number of initial infectives: I0 
% Note! I0 should be positive!
% number of initial recovered: R0
% number of initial diagnosed: D0
% trace not only those current infectious individuals but also those naturally recovered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a list of individuals who are APP users/ Non-App-Users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_App_users=[];
list_non_App_users=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decide whether individual using App or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
 for i=1:n
     
           rd_app = rand;
           
           if rd_app < p_A
               
            list_App_users = [list_App_users, i];
           
           else
               
               list_non_App_users =[list_non_App_users, i];
               
           end

     
 end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a table for the NUMBER of S,I,R,D at each time of event t: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of susceptibles at time t 
St=S0;
% the number of infectious individual at time t 
It=I0; 
% the number of naturally recovered individual at time t 
Rt=R0;
% the number of diagnosed individual at time t 
Dt=D0;
 
% number of individuals that have been infected during the epidemic
counts_inf=0;

% time array of events 
time=cell(1,1);
time{1,1}=0;
% array of the number of S at time
S_t=cell(1,1);
S_t{1,1}=St;
% array of the number of I at time
I_t=cell(1,1);
I_t{1,1}=It;
% array of the number of R at time
R_t=cell(1,1);
R_t{1,1}=Rt;
% array of the number of D at time
D_t=cell(1,1);
D_t{1,1}=Dt;


% Construct the table T_num:
T_number=table(time,S_t,I_t,R_t,D_t);



% construct a table of contacts and traced contacts:

% the column for labeled individuals:
individual=(1:n)';
% the list of contacts of each individual: 
contacts=cell(n,1);
% the list of manual reporting contacts of each individual: 
contacts_mct=cell(n,1);

T_contacts=table(individual,contacts,contacts_mct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a table for the LIST of S,I,R,D at each time of event t: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% randomly choose one of n individuals to be the index case

index=randperm(n,1);
index=index(1);
 
 
% list of Infectives 
% 
list_I=cell(1,1);
list_I{1,1}=[index];

%list of Susceptibles
% choose the other n-1 individuals except index to be susceptible
list_S=cell(1,1);
list_S{1,1}=[1:(index-1),(index+1):n];

% list of Recovered naturally
list_R=cell(1,1);

% list of Diagnosed
list_D=cell(1,1);


% Construct the table T_list:
T_list=table(time,list_S,list_I,list_R,list_D);


%% Start the Epidemic 

while It > 0 
    
% currently we have It infectious individuals 

% list of S,I,R,D from the last time step

old_list_S = list_S{size(list_S,1),1};
old_list_I= list_I{size(list_I,1),1};
old_list_R = list_R{size(list_R,1),1};
old_list_D = list_D{size(list_D,1),1};

% record rates of four events:

r_infection = beta *length(old_list_I)*(length(old_list_S)/n);


r_recovery = gamma*length(old_list_I);


r_diagnosis = delta*length(old_list_I);

% total rates:    
    r_total = r_infection+r_recovery+r_diagnosis;
    
% time until next event happens:

    dt=exprnd(1/r_total);




%% Now, randomly choose a uniform(0,1) number to decide which event exactly occurs
    
    rd=rand;

%% --------------------------------First, if new infection occurs--------------------------------------


if rd <= (r_infection/r_total) 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  NEW INFECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Randomly choose one of the remaining susceptibles to be infected:

       % random choose one position from the labeled susceptibles
       s_pos=randperm(length(old_list_S),1);

       % get exactly which susceptible it is
       infectee=old_list_S(s_pos);
      
% Randomly choose one of the infectives as the corr. infector 
      
       % random pick one of the infetives as the one who infected the ri1
       inf_pos=randperm(length(old_list_I),1);
 
       % get exactly which infective it is
        infector=old_list_I(inf_pos);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete the new infected from list of S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_list_S=old_list_S;

newinf_pos=find(new_list_S == infectee,1);
 
new_list_S(newinf_pos)=[];  
     
St=St-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the new infected to list of I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
new_list_I=[old_list_I,infectee];
 
% Record the infection contact between them
        
contacts{infectee,1}=[contacts{infectee,1},infector];
contacts{infector,1}=[contacts{infector,1},infectee];
        
It=It+1;     

counts_inf=counts_inf+1;

% Decide whether this contact would be manually traced:
        
           rdm = rand;
 
       
if rdm<p
       
            contacts_mct{infectee,1}=[infector,contacts_mct{infectee,1}];

            contacts_mct{infector,1}=[contacts_mct{infector,1},infectee];
else
               
end

% update the list of D
list_D{size(list_D,1)+1,1} = old_list_D;

% update the list of I
list_I{size(list_I,1)+1,1} = new_list_I;

% update the list of S
list_S{size(list_S,1)+1,1} = new_list_S;

% update the list of R 
list_R{size(list_R,1)+1,1} = old_list_R;

% record the number of I 

I_t{size(I_t,1)+1,1}=It;

% record the number of D
D_t{size(D_t,1)+1,1}=Dt;

% record the number of S

S_t{size(S_t,1)+1,1}=St;

% record the number of R
R_t{size(R_t,1)+1,1}=Rt;

%% ------------------------------------if Natural Recovery occurs ----------------------------------
    elseif rd > (r_infection/r_total) && rd <= ((r_infection+r_recovery)/r_total)
%  Randomly choose one of the remaining infectives to be naturally recovered:

      recover_pos=randperm(length(old_list_I),1);
      newrecover=old_list_I(recover_pos);
      Rt=Rt+1;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete the recovered from list of I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_list_I=old_list_I;
        
 
new_list_I(recover_pos)=[];  
     

It=It-1;       
 
% update the list of R

new_list_R=[old_list_R,newrecover];

list_R{size(list_R,1)+1,1} = new_list_R; 

% update the list of D
list_D{size(list_D,1)+1,1} = old_list_D;

% update the list of I
list_I{size(list_I,1)+1,1} = new_list_I;

% update the list of S
list_S{size(list_S,1)+1,1} = old_list_S;

% record the number of I 

I_t{size(I_t,1)+1,1}=It;

% record the number of D
D_t{size(D_t,1)+1,1}=Dt;

% record the number of S

S_t{size(S_t,1)+1,1}=St;

% record the number of R
R_t{size(R_t,1)+1,1}=Rt;


    
%% ------------------------------------ if Diagnosis occurs and tracing at same time----------------------------------     
else
 %%  Randomly choose one of the remaining infectives to be diagnosed: 
      
      Diag=[];    

   diagnose_pos=randperm(length(old_list_I),1);
    
   % record the number of D
   Dt=Dt+1;
   
   new_list_D=[old_list_D,old_list_I(diagnose_pos)];
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Delete the diagnosis from list of I
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new_list_I=old_list_I;
    new_list_I(diagnose_pos)=[];  
    It=It-1;       
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For those new diagnosed individuals:  CONTACT TRACING  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
% construct a list of individuals whose contacts will be traced.
Diag = [Diag, old_list_I(diagnose_pos)];

OldDiag=Diag;

% while the list is not empty:

while ~isempty(Diag) 
             
    
              if  It == 0 % if no infectives, stop the loop
                      break
              end
    
             NewDiag=[]; % ready to construct a new list of individuals whose contacts will be traced.
       
%------for each individual i from list Diag:
       
for i = 1:length(Diag)
           
           if  It == 0 % if no infectives, stop the loop
                      break
           end
         
          
%------- if this i-th individual is app-user

if ~isempty(find(list_App_users == Diag(i)))
               
%------- for each j-th contact of the  i-th individual :

            for j=1:length(contacts{Diag(i),1})
                
%-------- if the j-th contact is an app-user   

               con_j_i = contacts{Diag(i),1}(j);
               
               if  ~isempty(find(list_App_users == contacts{Diag(i),1}(j)))
                   
%-------- if the j-th contact has not been already diagnosed, then continue it with next follow steps: 

               if isempty(find(new_list_D == contacts{Diag(i),1}(j)))
                  
%-------- then there can be two cases: either it is still infectious or it is already naturally recovered: 
             
%-------- if the j-th contact is still infectious, then:

               if  ~isempty(find(new_list_I == contacts{Diag(i),1}(j)))
                  
                  %fprintf('\n the j contact is still infectious? %d\n',~isempty(find([labelI{:}] == symcon{Diag(i),1}(j))));
                 
                  % find out at which position of the list of I:
                  
                   i_j=find(new_list_I == contacts{Diag(i),1}(j));
                  
                  % Remove this j-th contact from infectives list of I:
                  new_list_I(i_j)=[];
                  It = It-1;
                  
                  % Add this j-th contact to the diagnosed list:
                  
                  new_list_D= [new_list_D,contacts{Diag(i),1}(j)];
                  
                  Dt = Dt+1;
                  
                  % Add this to be the new list NewDiag
                  
                  NewDiag=[NewDiag, contacts{Diag(i),1}(j)];
                  
                  if  It == 0 % if no infectives, stop the loop
                      break
                  end
                  
%-------- otherwise, this j-th contact is already naturally recovered, and we still tracing own contacts.
              
               else
                   
               %Since we still tracing its own contacts, add this to 
                  % the new list
                  
                  NewDiag=[NewDiag, contacts{Diag(i),1}(j)];
     
               end %-------- Row  

               else %-------- Row 
                   
               end %-------- Row 
               
               else
      
               end
               
           end %-------- Row
           
%-------- if the app-user Diag(i) can trigger manual tracing  

           if ~isempty(contacts_mct{Diag(i),1}) %-------- Diag(i) can trigger manual tracing
    
for j=1:length(contacts_mct{Diag(i),1})
%-------- if the j-th contact has not been already diagnosed, then continue it with next follow steps: 

if isempty(find(new_list_D == contacts_mct{Diag(i),1}(j)))   
%-------- then there can be two cases: either it is still infectious or it is already naturally recovered: 
             
%-------- if the j-th contact is still infectious, then:

if  ~isempty(find(new_list_I == contacts_mct{Diag(i),1}(j)))  
    
                  i_j=find(new_list_I == contacts_mct{Diag(i),1}(j));
                  
                  % Remove this j-th contact from infectives list of I:
                  new_list_I(i_j)=[];
                  It = It-1;
                  
                  % Add this j-th contact to the diagnosed list:
                  
                  new_list_D= [new_list_D,contacts_mct{Diag(i),1}(j)];
                  
                  Dt = Dt+1;
                  
                  % Add this to be the new list NewDiag
                  
                  NewDiag=[NewDiag, contacts_mct{Diag(i),1}(j)];
                  
                  if  It == 0 % if no infectives, stop the loop
                      break
                  end    
    
else %-------- otherwise, this j-th contact is already naturally recovered, and we still tracing own contacts.
                 NewDiag=[NewDiag, contacts_mct{Diag(i),1}(j)];
end

else
end

end                 

           else
           end
 
%-------- else Diag(i) is the non-app-user if it can trigger manual tracing

elseif ~isempty(contacts_mct{Diag(i),1}) 
    
for j=1:length(contacts_mct{Diag(i),1})
%-------- if the j-th contact has not been already diagnosed, then continue it with next follow steps: 

if isempty(find(new_list_D == contacts_mct{Diag(i),1}(j)))   
%-------- then there can be two cases: either it is still infectious or it is already naturally recovered: 
             
%-------- if the j-th contact is still infectious, then:

if  ~isempty(find(new_list_I == contacts_mct{Diag(i),1}(j)))  
    
                  i_j=find(new_list_I == contacts_mct{Diag(i),1}(j));
                  
                  % Remove this j-th contact from infectives list of I:
                  new_list_I(i_j)=[];
                  It = It-1;
                  
                  % Add this j-th contact to the diagnosed list:
                  
                  new_list_D= [new_list_D,contacts_mct{Diag(i),1}(j)];
                  
                  Dt = Dt+1;
                  
                  % Add this to be the new list NewDiag
                  
                  NewDiag=[NewDiag, contacts_mct{Diag(i),1}(j)];
                  
                  if  It == 0 % if no infectives, stop the loop
                      break
                  end    
    
else %-------- otherwise, this j-th contact is already naturally recovered, and we still tracing own contacts.
                 NewDiag=[NewDiag, contacts_mct{Diag(i),1}(j)];
end

else
end

end                 

                  
          
                  

end 


end      %--------  the tracing of all from the list Diag stops

 %-------- Now, we compare the new list NewDiag and the old list OldDiag:
       
       % if some elements in the NewDiag are also in the OldDiag, then
       % there is no need to trace them again from beginning, so delete
       % them from the NewDiag:
       
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


% update the list of D
list_D{size(list_D,1)+1,1} = new_list_D;

% update the list of I
list_I{size(list_I,1)+1,1} = new_list_I;

% update the list of S
list_S{size(list_S,1)+1,1} = old_list_S;

% update the list of R 
list_R{size(list_R,1)+1,1} = old_list_R;

% record the number of I 

I_t{size(I_t,1)+1,1}=It;

% record the number of D
D_t{size(D_t,1)+1,1}=Dt;

% record the number of S

S_t{size(S_t,1)+1,1}=St;

% record the number of R
R_t{size(R_t,1)+1,1}=Rt;
   
    
 
end

% record the time of event: 

% move to next time step 
time{size(time,1)+1,1}=time{size(time,1),1}+dt;



end %------------- end of while positive number of I 

% update the tables 
T_list=table(time,list_S,list_I,list_R,list_D);
save('T_list.mat', 'T_list')

T_contacts=table(individual,contacts);
save('T_contacts.mat', 'T_contacts')

T_number=table(time,S_t,I_t,R_t,D_t);
save('T_number.mat', 'T_number')

end %------------- end of function