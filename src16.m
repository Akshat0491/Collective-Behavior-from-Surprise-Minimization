%something still seems off!



%need to systematically figure out how parameters affect the simulation

%in current implementation, it looks at all other agents, but acts based on
%neightbours withing a particular radius

%can make many things into functions



tic

clc
clear
close all
%% Clearing memory
%% Parameters

dt=0.1;
k=0.01; %learning rate
Nitr=dt/k; %number of times belief is updated
T=50;    %in seconds
time=0:dt:T;
Nt=length(time); %Number of times
R=50; %what if this is not constant?
n=2;
Na=n; %Number of Agents
s=10;
S=20;
wanna_save=1;
have_noise=1;

theta0=-90; %degrees
max_theta=270;
L=6; %initialising sensory sectors
delta_theta=(max_theta-theta0)/L;
theta_max=theta0+(L-1)*delta_theta;
thetas=theta0:delta_theta:theta_max;

alpha1=.5;

eta=30;
            %noise parameters
         Lambda_w = 0.1;
         Gamma_w=1;
         Gamma_vecw= Gamma_w*ones(L,1);
         sigma_w= construct_sigma(Lambda_w,Gamma_vecw); %generating noise



          Lambda_z =1; %paramaters for sensory noise
          Gamma_z=1;
Gamma_vecz= Gamma_z*ones(L,1); %paramaters for noise
sigma_z= construct_sigma(Lambda_z,Gamma_vecz); %generating sensory noise

%%


pos = zeros(Na,2,Nt); %initialise r
vel = zeros(Na,2,Nt); %initialise rdot
acc = zeros(Na,2,Nt); %initialise rdotdot
jer = zeros(Na,2,Nt); %initialise rdotdotdot



%placing uniformly
d=s/(n-1);
for i=1:Na
    pos(i,:,1)=[-s/2 +  (mod(i-1,n))*d,-s/2 + (floor((i-1)/n))*d];
end
% pos(:,:,1) = s*rand(Na,2)-10; %at t=0, randomly placed in the square of sXs;  %how about placing them uniformly?
scatter(pos(:,1,1),pos(:,2,1),"filled")



mue_dists         = zeros(L,Na,Nt);
mue_dists_dot     = zeros(L,Na,Nt);
mue_dists_dot_dot = zeros(L,Na,Nt);
mue_dists_dot_dot_dot = zeros(L,Na,Nt);


Fl=zeros(Na,Nt); %initialising free energies




hidden_dists             = zeros(L,Na,Nt);
hidden_dists_dot         = zeros(L,Na,Nt); %initialised hidden environmental states
hidden_dists_dot_dot     = zeros(L,Na,Nt);
hidden_dists_dot_dot_dot = zeros(L,Na,Nt);





sense_dists         = zeros(L,Na,Nt);
sense_dists_dot     = zeros(L,Na,Nt); %initialised distances, can sense only distances
sense_dists_dot_dot = zeros(L,Na,Nt);


delta_rhat=zeros(L,2,Na,Nt); %for action update

%for now nois parameters remains same! can work on updating them later

%noise in dynamics



PI_z=inv(sigma_z);
PI_w=inv(sigma_w); 





angles=zeros(Na,Na,Nt);










% @t=1
% mue_pos(:,:,:,1)=eta*ones(L,2,Na); %belief of the agent, born with belief to be 2.5m apart 
% 



for i=1:Na

    angles(:,i,1)=atand((pos(:,2,1)-pos(i,2,1))./((pos(:,1,1)-pos(i,1,1))));
    
    for w=thetas
        agents=find(angles(:,i,1)>=w & angles(:,i,1)<=(w+delta_theta));
            
        if ~isempty(agents)
        % hidden_pos(thetas==w,:,i,1)=mean(pos(agents,:,1)-pos(i,:,1));
        % hidden_vel(thetas==w,:,i,1)=mean(vel(agents,:,1)-vel(i,:,1));
        % hidden_acc(thetas==w,:,i,1)=mean(acc(agents,:,1)-acc(i,:,1));
         w_tilde = mvnrnd(zeros(3*L,1),sigma_w)';
         % w_tilde_y = mvnrnd(zeros(3*L,1),sigma_w)';
         % w_tilde=[w_tilde_x, w_tilde_y];
if ~have_noise
    w_tilde=zeros(3*L,1);
end
          hidden_dists(thetas==w,i,1)=mean( get_dists(pos(agents,:,1)-pos(i,:,1)) );

          hidden_dists_dot(:,i,1)         = (- alpha1 * (hidden_dists(:,i,1)-eta))       + w_tilde(    1:  L);     
          hidden_dists_dot_dot(:,i,1)     = (- alpha1 * hidden_dists_dot(:,i,1))         + w_tilde(  L+1:2*L);           %initiated the derivatives of hidden states. right now, there's equal noise in x and y directions, can see how to change later
          hidden_dists_dot_dot_dot(:,i,1) = (- alpha1 * hidden_dists_dot_dot(:,i,1))     + w_tilde(2*L+1:3*L);
        end



    end


    z_tilde = mvnrnd(zeros(3*L,1),sigma_z)';


    sense_dists(:,i,1)         = hidden_dists(:,i,1) + z_tilde(    1:  L);  %g(x)=x; statement will have to be used repetitively at iterations, explore abstracting it into a fn object later
    sense_dists_dot(:,i,1)     = hidden_dists(:,i,1) + z_tilde(  L+1:2*L); %same as above
    sense_dists_dot_dot(:,i,1) = hidden_dists(:,i,1) + z_tilde(2*L+1:3*L);
    

    % mue_dists(:,i,1)         = sense_dists(:,i,1);
    % mue_dists_dot(:,i,1)     = sense_dists_dot(:,i,1); %keep on if this is how you want initial beliefs to be
    % mue_dists_dot_dot(:,i,1) = sense_dists_dot_dot(:,i,1);


    ez=[(sense_dists(:,i,1) - mue_dists(:,i,1));...
    (sense_dists_dot(:,i,1) - mue_dists_dot(:,i,1));...
    (sense_dists_dot_dot(:,i,1) - mue_dists_dot_dot(:,i,1))];

    ew=[ mue_dists_dot(:,i,1) + alpha1*(mue_dists(:,i,1) - eta); %be cautious about abs
         mue_dists_dot_dot(:,i,1) + alpha1*mue_dists_dot(:,i,1); 
         alpha1*mue_dists_dot_dot(:,i,1)];

    % Fl(i,1)= 0.5 * sum(sum((((ez')*(PI_z)*(ez)) + ((ew')*(PI_w)*(ew)) + (3*L*log(2*pi)))));
    Fl(i,1)= 0.5 * ((ez')*(PI_z)*(ez)) + ((ew')*(PI_w)*(ew)) + (3*L*log(2*pi));



    %action initialise yaha bhi krna hai

end

% mue_pos(:,:,:,1)=eta*ones(L,2,Na);


%%












for t=2:Nt %simultation loop
    
    ux=find((abs(pos(:,1,t-1))>S));
    uy=find((abs(pos(:,2,t-1))>S));

    vel(:,:,t) = vel(:,:,t-1) + dt * acc(:,:,t-1);
    if ~isempty(ux)
        vel(ux,1,t)=-1*vel(ux,1,t); %reflection about the boundary
    end
    if ~isempty(uy)
        vel(uy,2,t)=-1*vel(uy,2,t); %reflection about the boundary
    end
    


    pos(:,:,t) = pos(:,:,t-1) + 10*dt * (vel(:,:,t)./get_dists(vel(:,:,t))); %same ) speeds for all, is it justified!? time evolution of position
    

    



    for i=randperm(Na)
      % eta=repmat([2.5 2.5],L,1); %currently same for all, can make it diff for each agent later



         w_tilde = mvnrnd(zeros(3*L,1),sigma_w)';


if ~have_noise
    w_tilde=zeros(3*L,1);
end

          
      angles(:,i,t)=atand((pos(:,2,t)-pos(i,2,t))./((pos(:,1,t)-pos(i,1,t))));
    
    for w=thetas
        agents=find(angles(:,i,t)>=w & angles(:,i,t)<=(w+delta_theta));
            
        if ~isempty(agents)
         w_tilde = mvnrnd(zeros(3*L,1),sigma_w)';

if ~have_noise
    w_tilde=zeros(3*L,1);
end
          hidden_dists(thetas==w,i,t)=mean( get_dists(pos(agents,:,t)-pos(i,:,t)) );

          hidden_dists_dot(:,i,t)         = (- alpha1 * (hidden_dists(:,i,t)-eta))       + w_tilde(    1:  L);     
          hidden_dists_dot_dot(:,i,t)     = (- alpha1 * hidden_dists_dot(:,i,t))         + w_tilde(  L+1:2*L);           %initiated the derivatives of hidden states. right now, there's equal noise in x and y directions, can see how to change later
          hidden_dists_dot_dot_dot(:,i,t) = (- alpha1 * hidden_dists_dot_dot(:,i,t))     + w_tilde(2*L+1:3*L);
        end



    end






%% sensory states
%An agent senses the position of all the other agents for now
% y = x + w, i.e g(x)=x, refernce rom the paper






    z_tilde = mvnrnd(zeros(3*L,1),sigma_z)';
if ~have_noise
    z_tilde=zeros(3*L,1);
end

    sense_dists(:,i,t)         = hidden_dists(:,i,t) + z_tilde(    1:  L);  %g(x)=x; statement will have to be used repetitively at iterations, explore abstracting it into a fn object later
    sense_dists_dot(:,i,t)     = hidden_dists(:,i,t) + z_tilde(  L+1:2*L); %same as above
    sense_dists_dot_dot(:,i,t) = hidden_dists(:,i,t) + z_tilde(2*L+1:3*L);






        Dmue=[mue_dists_dot(:,i,t-1);mue_dists_dot_dot(:,i,t-1);mue_dists_dot_dot_dot(:,i,t-1)];

                


        gez=[eye(L), zeros(L), zeros(L);
        zeros(L),eye(L),zeros(L);
        zeros(L),zeros(L),eye(L)];


    ez=[sense_dists(:,i,t) - mue_dists(:,i,t-1);...
        sense_dists_dot(:,i,t) - mue_dists_dot(:,i,t-1);...
        sense_dists_dot_dot(:,i,t) - mue_dists_dot_dot(:,i,t-1)];

    gez=gez' * PI_z * ez;


    few=[-alpha1*eye(L), zeros(L), zeros(L);
        zeros(L),-alpha1*eye(L),zeros(L);
        zeros(L),zeros(L),-alpha1*eye(L)];

    ew=[ mue_dists_dot(:,i,t-1) + alpha1*(mue_dists(:,i,t-1) - eta); %be cautious about abs
        mue_dists_dot_dot(:,i,t-1) + alpha1*mue_dists_dot(:,i,t-1); 
        alpha1*mue_dists_dot_dot(:,i,t-1)];


    few=few' * PI_w * ew;

    dtew = - [zeros(L),zeros(L),zeros(L);
              eye(L),zeros(L),zeros(L);
              zeros(L),eye(L),zeros(L)]*ew;


              for w=thetas
                        agents=find(angles(:,i,t)>=w & angles(:,i,t)<=(w+delta_theta));
                        if ~isempty(agents)
          Nin=find( get_dists(pos(agents,t)-pos(i,t))<R & get_dists(pos(agents,t)-pos(i,t))>0); %interaction set, but wait, what if it is empty??.....%optimise it later
          if ~isempty(Nin) %for now, interact only if something is nearby
          K=length(Nin);

          delta_rhat(w==thetas,:,i,t)= mean((pos(agents(Nin),:,t)-pos(i,:,t)) ./ get_dists((pos(agents(Nin),:,t)-pos(i,:,t)))); %average unit vectors for each agent in interaction set

          end
                        end

          end

              if t> 0 %(20/100)*Nt %action updates comes into play after 1/5th time
          % acc(i,:,t)=2*Gamma_z*Lambda_z*Lambda_z*sum((sense_vel(:,:,i,t)-mue_dists_dot(:,:,i,t)).*delta_rhat(:,:,i,t)); %might go wrong!
          acc(i,1,t)=2*Gamma_z*Lambda_z*Lambda_z*(   ((sense_dists(:,i,t)-mue_dists(:,i,t-1))' ) * delta_rhat(:,1,i,t));
          acc(i,2,t)=2*Gamma_z*Lambda_z*Lambda_z*(   ((sense_dists(:,i,t)-mue_dists(:,i,t-1))' ) * delta_rhat(:,2,i,t));
              end    %action update before belief update


    mue_dists(:,i,t)         = mue_dists(:,i,t-1)        + k*(Dmue(1:L)       + gez(1:L)       + few(1:L)       + dtew(1:L));
    mue_dists_dot(:,i,t)     = mue_dists_dot(:,i,t-1)    + k*(Dmue(L+1:2*L)   + gez(L+1:2*L)   + few(L+1:2*L)   + dtew(L+1:2*L));
    mue_dists_dot_dot(:,i,t) = mue_dists_dot_dot(:,i,t-1)+ k*(Dmue(2*L+1:3*L) + gez(2*L+1:3*L) + few(2*L+1:3*L) + dtew(2*L+1:3*L));

    % for tm=1:Nitr-1 %Nitr
    %         mue_pos(:,:,i,t) = mue_pos(:,:,i,t)+ k*(Dmue(1:L,:)       + gez(1:L,:)       + few(1:L,:)       + dtew(1:L,:));
    % mue_vel(:,:,i,t) = mue_vel(:,:,i,t)+ k*(Dmue(L+1:2*L,:)   + gez(L+1:2*L,:)   + few(L+1:2*L,:)   + dtew(L+1:2*L,:));
    % mue_acc(:,:,i,t) = mue_acc(:,:,i,t)+ k*(Dmue(2*L+1:3*L,:) + gez(2*L+1:3*L,:) + few(2*L+1:3*L,:) + dtew(2*L+1:3*L,:));
    % end
    % 
          



    Fl(i,t)= 0.5 * (((ez')*(PI_z)*(ez)) + ((ew')*(PI_w)*(ew)) + (3*L*log(2*pi)));
    
    end






end

% plot(time,Fl(1,:))
% plot(time,squeeze(pos(1,1,:)))










    % acc(:,:,t+1)=acc(:,:,t) + dt*(jer(:,:,t+1));
    % vel(:,:,t+1)=vel(:,:,t) + dt*(acc(:,:,t+1)); 
    % pos(:,:,t+1)=pos(:,:,t) + dt*(vel(:,:,t+1));


    
%     pos_t=pos(:,:,t);
%     vel_t=vel(:,:,t);
%     acc_t=acc(:,:,t);
%     jer_t=jer(:,:,t);
% 
% 
%     for i=1:Na %optimization chance: replace with parrallel for loop iterating for every individual
%           rel_pos=pos_t-pos(i,:,t); %hidden states x, as per the paper
%           L=Na
%           alpha1=1;
%           eta=repmat([2.5 2.5],L,1); %currently same for all, can make it diff for each agent later
% 
% 
%          Lambda_w = 1;
%          Gamma_vecw= ones(L,1);
%          sigma_w= construct_sigma(Lambda_w,Gamma_vecw); %generating noise
%          w_tilde = mvnrnd(zeros(3*L,1),sigma_w)';
% 
% 
% 
%           rel_vel = (- alpha1 * (rel_pos-eta)) + w_tilde(1:L);     
%           rel_acc = (- alpha1 * rel_vel) + w_tilde(L+1:2*L);           %initiated the derivatives of hidden states. right now, there's equal noise in x and y directions, can see how to change later
%           rel_jer = (- alpha1 * rel_acc) + w_tilde(2*L+1:3*L);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% sensing
% 
% 
% Lambda_z = 1; %paramaters for noise
% Gamma_vecz= ones(Na,1); %paramaters for noise
% sigma_z= construct_sigma(Lambda_z,Gamma_vecz); %generating sensory noise
% z_tilde = mvnrnd(zeros(3*Na,1),sigma_z)';
% 
% 
% sense_pos = rel_pos + z_tilde(1:Na);  %g(x)=x; statement will have to be used repetitively at iterations, explore abstracting it into a fn object later
% sense_vel = rel_vel + z_tilde(Na+1:2*Na); %same as above
% sense_acc = rel_acc + z_tilde(2*Na+1:3*Na);
% 
% 
% 
% 
% 
%     end
% 
% 

% end

% scatter(time,Fl(2,:),'g')

% scatter(time,Fl(1,:))
% scatter(time,Fl(3,:))
% scatter(time,Fl(4,:))
% plot(squeeze(mue_pos(1,1,2,:)),Fl(2,:))
% plot(squeeze(hidden_pos(1,1,2,:)),Fl(2,:))

% plot(time,reshape(mue_pos(1,1,2,:),[1 Nt]))
toc1=toc
% plot(time,reshape(mue_pos(2,1,1,:),[1 Nt]),'w')
% end

% function get_neighours




%% Plotting and saving

    %took from chat gpt
% Set up video writer

figure(1)
plot(time,Fl([1,floor(Na/2),Na],:))
title("Free Energy")
xlabel("Time")
ylabel("Free Energy")


figure(2)
plot(time,squeeze(sqrt(sum((pos(2,:,:)-pos(1,:,:)).^2,2))))

title("Inter Agent Distance")
xlabel("Time")
ylabel("Distance")


if wanna_save
output_dir ="D:\Projects\Summer 2025\Surprise Minimisation\outputs"; % replace with your desired path

now_dt = datetime('now');
now_dt.TimeZone = 'local';
now=string(now_dt, 'd-MMM-y HH:mm:ss');
now=split(now,":");
now=join(now,"_");
filename = fullfile(output_dir, sprintf('Agent Simulation_%s', ...
now));


outputVideo = VideoWriter(filename,'MPEG-4');
outputVideo.FrameRate = 15; % Adjust as needed
% outputVideo.FileFormat='mp4';
open(outputVideo);

fig = figure('Color','k'); % black background for visual clarity

for t = 1:Nt-10
clf
%scatter(pos(:,1,t), pos(:,2,t), 10, 'w', 'filled'); % white dots
quiver(pos(:,1,t), pos(:,2,t), 2*vel(:,1,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps), 2*vel(:,2,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps),"off","Marker",".","ShowArrowHead","on")
% xlabel(sprintf('X — sensed: %.2f, belief: %.2f', ...
%     mean(sense_dists(:,:,t), 'all'), ...
%     mean(mue_dists(:,:,t), 'all')));

ylabel("Y");
xlim([-S S]);
ylim([-S S]);

if have_noise
title(sprintf('Time %.2f sec,N:%.f,\n dt:%.2f, k:%.2f, R:%.2f, vision:full,\n alpha:%f, eta:%.2f, Lambda_w:%.2f,Lambda_z:%.2f, Gamma_w:%.2f, Gamma_z: %.2f', (t-1)*dt,Na,dt,k,R,alpha1,eta,Lambda_w,Lambda_z,Gamma_w,Gamma_z), 'Color', 'w');
else
title(sprintf('Without noise Time %.2f sec,N:%.f,\n dt:%.2f, k:%.2f, R:%.2f, vision:full,\n alpha:%f, eta:%.2f, Lambda_w:%.2f,Lambda_z:%.2f, Gamma_w:%.2f, Gamma_z: %.2f', (t-1)*dt,Na,dt,k,R,alpha1,eta,Lambda_w,Lambda_z,Gamma_w,Gamma_z), 'Color', 'w');
end
set(gca, 'Color', 'k'); % black background for axes
drawnow;

% Capture the plot as a frame
frame = getframe(fig);
writeVideo(outputVideo, frame);
end    
    close(outputVideo);
    disp('Video saved successfully.');
    filename
    end


 



function Sigma = construct_sigma(lambda, Gamma_vec)
sigma_s=[1,0,-0.5/(lambda^2);
         0,0.5/(lambda^2), 0;
         -0.5/(lambda^2),0,0.75/(lambda^4)];
sigma_t=diag(Gamma_vec);
Sigma=kron(sigma_t,sigma_s);
end


function  dists = get_dists(pos)
dists = sqrt(sum(pos.^2,2))+eps;
end