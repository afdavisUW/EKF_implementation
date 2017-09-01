function [sim] = smoother_RTS_Lin( func, y_meas )
% Function to execute the EKF and RTS smoother
%   Certain sections can be commented out in order to allow just the EKF to
%   run and not the entire RTS

global Param Time Data Test %Func


tic
% extracting dynamics and output functions
% f = func.fODE; % used because it calculates correct particle veloctity
h = func.h;
F = func.F;
H = func.H;

y.meas = y_meas'; % smoother function is expecting a row vector not a column vector

%% First the EKF must run
% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initial values
x.ini = Param.LinIC.est;
P.ini = Param.LinIC.P;
% vector sizes
x.apri = zeros(size(x.ini));
x.apri = x.ini
x.apost = zeros(size(x.ini));
P.apri = P.ini;

% null initializations
Time.hat = [];
x.hat = [];
P.vals = [];
Param.F = [];

% Begin EKF Loop
for k = 1:length(Time.t)-1

if k > 100 %573
    bob = 1;
end

    if mod(k, 10) == 0
        display(['Progress: ',num2str(k),'/',num2str(length(Time.t)-1),': C_D = ',num2str(x.apri(3))])
    end
    
%     Gain
    K = P.apri*H(x.apri)'/(H(x.apri)*P.apri*H(x.apri)'+Param.R);
    Param.savedK(:,:,k) = K;
    
%     Update
    if k == 1 % to account for the initial xhat
        x.apost = x.apri;
        P.apost = P.apri;
    else
        x.apost = x.apri' + (K*(y.meas(:,k) - h(x.apri')))';
%         P.apost = (eye(Param.numStates) - K*H(x.apri))*P.apri; % standard covariance update
%         PartialIden = zeros(Param.numStates); PartialIden(1,1) = 1; PartialIden(2,2) = 1; PartialIden(3,3) = 1;
%         P.apost = (PartialIden - K*H(x.apri))*P.apri; % standard covariance update

%         Joseph form update
        P.apost = (eye(Param.numStates) - K*H(x.apri))*P.apri...
            *(eye(Param.numStates)-K*H(x.apri))' + ...
            K*Param.R*K';
%         P.apost = (PartialIden - K*H(x.apri))*P.apri...
%             *(PartialIden-K*H(x.apri))' + ...
%             K*Param.R*K';
        Ptemp = P.apost(1:3,1:3);
        P.apost = zeros(Param.numStates);
        P.apost(1:3,1:3) = Ptemp;

    end

%     Propigation
    Time.sim = [Time.t(k):Time.dtSim:Time.t(k+1)-Time.dtSim];
    
%     Global test variable to duplicate breaking
Test.time = Time.sim;
Test.x = x.apost;

    
    [~,Xsim] = ode45(func.fODE, Time.sim, x.apost,Param.options);
    Time.hat = [Time.hat Time.sim];
    x.hat = [x.hat; Xsim];
    
    for j = 1:length(Xsim)
        [U.posTemp,U.velTemp,U.accTemp] = particleDynamics_spectrum(Time.sim(j));
        Param.Fsim(:,:,j) = F(Xsim(j,:), U.posTemp, U.velTemp, U.accTemp);
    end
        Param.F = cat(3, Param.F, Param.Fsim);

        
%     Using ode solver to solve for the covariance
    [~,PsimCol] = ode45('covariance_EKF', Time.sim, P.apost(:),Param.options);
%     [~,PsimCol] = ode45('CovarianceSolver', Time.sim, P.apost(:));

    P.vals = [P.vals; PsimCol];
    
%     Initialize next step
    P.apri = reshape(PsimCol(end,:),Param.numStates,Param.numStates)';
    x.apri = x.hat(end,:)';
    

end % end EKF for loop

%     Transpose to maintain consistency with reshapes as formattted below
    x.hat = x.hat';
    P.vals = P.vals';
    
%     Compute Output
for k = 1:length(x.hat)
    y.hat(:,k) = h(x.hat(:,k)');
end

%     Compute 3 Sigma Error Bounds
P.vals = reshape(P.vals, Param.numStates, Param.numStates, length(P.vals));

variance = zeros(Param.numStates, length(P.vals));
sig3 = zeros(Param.numStates, length(P.vals));
for k = 1:Param.numStates
    variance(k,:) = P.vals(k,k,:);
    sig3(k,:) = 3*variance(k,:).^(1/2);
end

% Save EKF data
EKF.time = Time;
EKF.x = x;
EKF.y = y;
EKF.P = P;
EKF.sig3 = sig3;

toc
%% Second Loop is the RTS smoother
if Param.RTSBool == 1 % excecute if RTS is active
    clear x P sig3 Param.Fsim

    % initialize smoother
    x.fin = EKF.x.hat(:,end);
    P.fin = EKF.P.vals(:,:,end);
    Param.x.f = timeseries(EKF.x.hat,EKF.time.hat,'NAME','xhat_EKF');

    % calculating the smoother gain
    for k = 1:length(EKF.P.vals)
        Param.K(:,:,k) = Param.G*Param.Q*Param.G'/EKF.P.vals(:,:,k);
    end

    % integrating the dynamics backwards in time for the smoother
    [~, xhat] = ode45(@(t, x) xdot_RTS(t, x, func.f), fliplr(Time.hat), x.fin); % simulating my xhat values
    x.hat = fliplr(xhat');
    toc
    % Integrate the covariance values
    [~,PsimCol] = ode45('covariance_RTS', fliplr(Time.hat), P.fin(:) ); % (:) colon notation reorients a matrix as a column vector
    % reorient then flip P.vals to be forward in time
    PsimCol = fliplr(PsimCol');
    P.vals = reshape(PsimCol, Param.numStates, Param.numStates, length(PsimCol));
    toc
    variance = zeros(Param.numStates, length(P.vals));
    sig3 = zeros(Param.numStates, length(P.vals));
    for k = 1:Param.numStates
        variance(k,:) = P.vals(k,k,:);
        sig3(k,:) = 3*variance(k,:).^(1/2);
    end

    % saving the RTS data
    RTS.time = Time;
    RTS.x = x;
    RTS.y = y;
    RTS.P = P;
    RTS.sig3 = sig3;

end

%% Save Simulation Data
sim.time = Time;
sim.EKF = EKF;
if Param.RTSBool == 1 % excecute if RTS is active
    sim.RTS = RTS;
end 
%% Remove extra parameter fields
fields = {'Fsim','savedK'};
Param = rmfield(Param,fields);



end

