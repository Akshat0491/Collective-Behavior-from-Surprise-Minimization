function ret=sense(Y_tilde,R_tilde,G,noise_params,t)
    %Y_tilde=cell(o,1)
    %Y_tilde{1}=cell(N,1)
    %Y_tilde{1}{1}=g(R_tilde{1})
    %decide whether to have noise or not on case by case basis
    % R_tilde : Nx2xT matrix, where N is the number of agents
    % G : Nx2 --> L function handle that takes R_tilde and returns the observation Y.
    % noise_params : parameters for the noise, if any, can be empty
    % Y_tilde : Observable states created by G(R)
    % t : current time step.
    o=length(Y_tilde);
    N=size(R_tilde{1},1);

    for p=1:o
        Y_tilde_p=Y_tilde{p};
        R_tilde_p=R_tilde{p}; 
        parfor i=1:N
            Y_tilde_p{i}(:,t)=G(R_tilde_p(:,:,t),i); %+noise(noise_params) %Y_tilde{p}{i} is the observation of agent i at time t, given the state of the world R_tilde{p}
        end
        Y_tilde{p}=Y_tilde_p; %update the cell array with the new observations
    end
    ret=Y_tilde; %return the updated observations

end


%write a function that takes R_tilde (the 2D Vector) and the Function G, noise parameters and returns the observation Y., 







