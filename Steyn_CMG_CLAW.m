
initW = [0.01;-0.01;0.005];
initAtt = [0.1; 0.2; 0.3];
initGama = deg2rad([45;-45;45;-45]);
% initGama = deg2rad([0;0;0;0]);
%[-60;60;120;-120]); % this should be scalar

gama_d = deg2rad([0;0;0;0]); % this should be scalar
initGama_dot = [0;0;0;0];  % this should be scalar
initOmega_wheels = [14.4;14.4;14.4;14.4];

initTheta = deg2rad(54.75);
N = 4;
gs_0(1:3,1) = [0;1;0];
gs_0(1:3,2) = [0;-1;0];
gs_0(1:3,3) = [1;0;0];
gs_0(1:3,4) = [-1;0;0];
Gs_0 = gs_0;

gg(1:3,1) = [cos(initTheta);0;sin(initTheta)];
gg(1:3,2) = [-cos(initTheta);0;sin(initTheta)];
gg(1:3,3) = [0; cos(initTheta);sin(initTheta)];
gg(1:3,4) = [0; -cos(initTheta);sin(initTheta)];
Gg = gg;
% gg(1:3,1) = [0;0;0];
% gg(1:3,2) = [0;0;0];
% gg(1:3,3) = [0;0;0];
% gg(1:3,4) = [0;0;0];
% Gg = gg;

gt_0(1:3,1) = cross(gg(1:3,1),gs_0(1:3,1));
gt_0(1:3,2) = cross(gg(1:3,2),gs_0(1:3,2));
gt_0(1:3,3) = cross(gg(1:3,3),gs_0(1:3,3));
gt_0(1:3,4) = cross(gg(1:3,4),gs_0(1:3,4));
Gt_0 = gt_0;

G_0 = [Gs_0;Gt_0;Gg];

tMax = 300;
delT = 0.1;
t = 0:delT:tMax;

I_sc = [86,0,0;0,85,0;0,0,113];
I_w_s = 0.1;
h_mag = I_w_s*100*[1;1;1];
wheelSpeedCurr = 100*[1;1;1;1];
L_ext = [0;0;0];
% L_ext = [0.1;0.1;-0.1];
u_s = [0;0;0;0];
u_g = [0;0;0;0];
hs = zeros(N,length(t));
Js_g = 0.13;
Jt_g = 0.02;
I_w_t = Jt_g;
Jg_g = 0.03;
J_g = diag([I_w_s,I_w_t,I_w_t]);  
beta = deg2rad(54.75);

% Td = 120;
P = 40; %(2*I(1,1))/Td;
K = 20; %(P(1,1)^2)/I(2,2);
K_gamaDot = [10; 10;10;10]; % based on p and k eqn % 3.82 (1/3 of P) based on 15s convergence time? Td = 15

rho = 0.8;
alpha = 0.001;

X = zeros(14,length(t));
X_dot = zeros(14,length(t));
X(:,1) = [initAtt;initW;initGama;initGama_dot];

% u = zeros(8,length(t));
omegaErr = zeros(3,length(t));
sigmaErr = zeros(3,length(t));

D1 = zeros(3,N);
D2 = zeros(3,N);
D3 = zeros(3,N);
D4 = zeros(3,N);

f = 0.03;
counter = 0;

for i = 1:(length(t)-1)
    
    % calculate sigma_R/N
    R_sigma_RN = (1/4)*[0.1*sin(f*t(i)); 0.2*cos(f*t(i));-0.3*sin(2*f*t(i))];
    
    % calculate sigma_R/N_dot by differentiating sigma_R/N
    R_sigma_RN_dot = (1/4)*[0.1*cos(f*t(i))*f; -0.2*sin(f*t(i))*f;-0.3*cos(2*f*t(i))*2*f];
    
    % calculate omega_R/N from Bmat of sigma_R/N and sigma_R/N_dot
    R_omega_RN = 4*(BmatMRP(R_sigma_RN)\R_sigma_RN_dot);
    
    % calculate sigma_R/N for next time step dt 
    R_sigma_RN_dt = (1/4)*[0.1*sin(f*t(i+1)); 0.2*cos(f*t(i+1));-0.3*sin(2*f*t(i+1))];
    
    % calculate sigma_R/N for next time step dt 
    R_sigma_RN_dot_dt = (1/4)*[0.1*cos(f*t(i+1))*f; -0.2*sin(f*t(i+1))*f;-0.3*cos(2*f*t(i+1))*2*f];
    
    R_omega_RN_dt = 4*(BmatMRP(R_sigma_RN_dt)\R_sigma_RN_dot_dt);
    
    R_omega_RN_dot = (R_omega_RN_dt- R_omega_RN)/delT;

    sigmaErr(:,i) = calcSigmaErr(X(1:3,i), R_sigma_RN);
    
    if norm(sigmaErr(:,i)) > 1
        sigmaErr(:,i) = -(sigmaErr(:,i)/(norm(sigmaErr(:,i))^2));
    end

    BR = MRP2C(sigmaErr(:,i));
    omega_RN_dot = BR*R_omega_RN_dot;
    omega_RN = BR*R_omega_RN;

    omegaErr(:,i) = (X(4:6,i)-omega_RN);
    omega_tilde_BN = [0, -X(6,i), X(5,i); ...
                  X(6,i), 0, -X(4,i); ...
                      -X(5,i), X(4,i), 0];
    
    gs_0 = G_0(1:3,1:4);
    gt_0 = G_0(4:6,1:4);
    gg   = G_0(7:9,1:4);
    gama = X(7:10,i);
    gamaDotCurr = X(11:14,i);
%     wheelSpeedCurr = X(15:18);
    
    h = h_mag(1)*[-cos(beta)*sin(gama(1))-cos(gama(2))+cos(beta)*sin(gama(3))-cos(gama(4)); ...
         -cos(gama(1))-cos(beta)*sin(gama(2))+cos(gama(3))+cos(beta)*sin(gama(4)); ...
         sin(beta)*sin(gama(1))+sin(beta)*sin(gama(2))+sin(beta)*sin(gama(3))+sin(beta)*sin(gama(4))];
%      
    A = h_mag(1)*[-cos(beta)*cos(gama(1)),    -sin(gama(2)),             cos(beta)*cos(gama(3)), -sin(gama(4)); ...
         -sin(gama(1)),              -cos(beta)*cos(gama(2)),    sin(gama(3)),           cos(beta)*cos(gama(4)); ...
         sin(beta)*cos(gama(1)),     sin(beta)*cos(gama(2)),    sin(beta)*cos(gama(3)),  sin(beta)*cos(gama(4))];
    h_dot = A*gamaDotCurr;
    
    for k = 1:N
        delAngle(k) = gama(k) - initGama(k); % subtracting from initial gama instead of previous gama because error builds up using previous gama instead of the absolute or constant value of initial gama. also this euation wouldnt work for previous gama. 
        gs(:,k) = [cos(delAngle(k))*gs_0(1:3,k) + sin(delAngle(k))*gt_0(1:3,k)];
        gt(:,k) = [-sin(delAngle(k)).*gs_0(1:3,k) + cos(delAngle(k)).*gt_0(1:3,k)];
    end
    
    omega_bn = X(4:6,i);
    omega_s_1 = gs(1:3,1)'*omega_bn;
    omega_t_1 = gt(1:3,1)'*omega_bn;
    omega_g_1 = gg(1:3,1)'*omega_bn;

    omega_s_2 = gs(1:3,2)'*omega_bn;
    omega_t_2 = gt(1:3,2)'*omega_bn;
    omega_g_2 = gg(1:3,2)'*omega_bn;

    omega_s_3 = gs(1:3,3)'*omega_bn;
    omega_t_3 = gt(1:3,3)'*omega_bn;
    omega_g_3 = gg(1:3,3)'*omega_bn;

    omega_s_4 = gs(1:3,4)'*omega_bn;
    omega_t_4 = gt(1:3,4)'*omega_bn;
    omega_g_4 = gg(1:3,4)'*omega_bn;

    omega_s = [omega_s_1; omega_s_2; omega_s_3;omega_s_4];
    omega_t = [omega_t_1; omega_t_2; omega_t_3;omega_t_4];
    omega_g = [omega_g_1; omega_g_2; omega_g_3;omega_g_4];
    BG_1 = [gs(1:3,1),gt(1:3,1),gg(1:3,1)];
    GB_1 = BG_1';
    BG_2 = [gs(1:3,2),gt(1:3,2),gg(1:3,2)];
    GB_2  = BG_2';
    BG_3 = [gs(1:3,3),gt(1:3,3),gg(1:3,3)];
    GB_3 = BG_3';
    BG_4 = [gs(1:3,4),gt(1:3,4),gg(1:3,4)];
    GB_4 = BG_4';

    J_b_1 = [BG_1]*[J_g]*[GB_1];
    J_b_2 = [BG_2]*[J_g]*[GB_2];
    J_b_3 = [BG_3]*[J_g]*[GB_3];
    J_b_4 = [BG_4]*[J_g]*[GB_4];
    J_b = [J_b_1,J_b_2,J_b_3,J_b_4];
    I = I_sc + J_b_1 + J_b_2 + J_b_3 + J_b_4;

    D0 = gs*I_w_s;
    h_nom = Js_g*initOmega_wheels;
    mu = 1e-9;
    u_cmg_part = 0;
% 
    for w = 1:N
%         D1(:,w) = ((I_w_s*wheelSpeedCurr(w)' + (Js_g/2)*omega_s(w))*gt(1:3,w)) + ((Js_g/2)*omega_t(w)*gs(1:3,w));
%         D2(:,w) = (0.5*Jt_g*(omega_t(w)*gs(1:3,w) + omega_s(w)*gt(1:3,w)));
%         D3(:,w) = (Jg_g*(omega_t(w)*gs(:,w) - omega_s(w)*gt(:,w)));
%         D4(:,w) = 0.5*(Js_g - Jt_g)*(gs(:,w)*gt(:,w)'*omega_RN + gt(:,w)*gs(:,w)'*omega_RN);
        
        D1(:,w) = ((I_w_s*wheelSpeedCurr(w)' + (I_w_s/2)*omega_s(w))*gt(1:3,w)) + ((I_w_s/2)*omega_t(w)*gs(1:3,w));
        D2(:,w) = (0.5*I_w_t*(omega_t(w)*gs(1:3,w) + omega_s(w)*gt(1:3,w)));
        D3(:,w) = (I_w_t*(omega_t(w)*gs(:,w) - omega_s(w)*gt(:,w)));
        D4(:,w) = 0.5*(I_w_s - I_w_t)*(gs(:,w)*gt(:,w)'*omega_RN + gt(:,w)*gs(:,w)'*omega_RN);
        
%         delta(i,w) = det((1/(h_nom(w))^2) * D1(1:3,w)*D1(1:3,w)');
        u_cmg_part = u_cmg_part + I_w_s*((wheelSpeedCurr(w)*omega_g(w)*gt(1:3,w)) - (wheelSpeedCurr(w)*omega_t(w)*gg(1:3,w))); 
    end
    D = D1 - D2 + D3 + D4; %full D matrix
%     D = D1 + D4; % [tracking problem: D matrix w/o insignificant terms like Jt, Jg, w_g, w_t, w_s
%     D = D1; % [regulation problem: D matrix w/o insignificant terms like Jt, Jg, w_g, w_t, w_s, omega_R/N

    %     Q = [D0, D];
%     W_s = (200*exp(-mu*delta))';
%     W_g = [1;1;1;1];
%     W = [            W_s(1), zeros(1,7); ...
%          zeros(1,1), W_s(2), zeros(1,6); ...
%          zeros(1,2), W_s(3), zeros(1,5); ...
%          zeros(1,3), W_s(4), zeros(1,4); ...
%          zeros(1,4), W_g(1), zeros(1,3); ...
%          zeros(1,5), W_g(2), zeros(1,2); ...
%          zeros(1,6), W_g(3), zeros(1,1); ...
%          zeros(1,7), W_g(4)];
   
%     u(:,i) = (W*Q'*inv(Q*W*Q'))*(-(-K*sigmaErr(:,i) - P*omegaErr(:,i) + I*(omega_RN_dot - omega_tilde_BN*omega_RN) + omega_tilde_BN*I*X(4:6,i)) +  u_cmg_part);    
%     u(:,i) = *(-(-K*sigmaErr(:,i) - P*omegaErr(:,i) + I_rw*(omega_RN_dot - (omega_tilde_BN*omega_RN)) + omega_tilde_BN*(I_rw*X(4:6,i) + Gs_0*hs(:,i))));      
    
% penrose method, pseudo ineverse steerign logic  and vadali control
%     gama_diff = gama_d - gama;
%     vadali_control = rho*(eye(4,4)-(A')*(inv(A*A' + alpha*eye(3,3))*A))*gama_diff;
%     u(:,i) = (A'*inv(A*A'))*-(-K*sigmaErr(:,i) - P*omegaErr(:,i) + I_sc*(omega_RN_dot - (omega_tilde_BN*omega_RN)) + omega_tilde_BN*I_sc*X(4:6,i) + omega_tilde_BN*h); %+ vadali_control;
% nakamura and hanafusa 
%     Lr = -K*sigmaErr(:,i) - P*omegaErr(:,i) + I*(omega_RN_dot - omega_tilde_BN*omega_RN) + omega_tilde_BN*I*X(4:6,i) +  u_cmg_part;

% [CLAW] nakamura and hanafusa
% u(:,i) = D'*inv(D*D' + alpha*eye(3,3))*(-Lr);

% [CLAW] penrose moore
% u(:,i) = D'*inv(D*D')*(-Lr);
%     Lr = -K*sigmaErr(:,i) - P*omegaErr(:,i) + I*(omega_RN_dot - omega_tilde_BN*omega_RN) + omega_tilde_BN*I*X(4:6,i) + omega_tilde_BN*h;
    Lr(:,i) = -K*sigmaErr(:,i) - P*omegaErr(:,i) + omega_tilde_BN*h;

    u(:,i) = A'*inv(A*A')*(-Lr(:,i)); % goes wild, bad claw and dyn behavior
%     using A matrix, maybe because not meant to track, only reg problem?

% [CLAW] vadali null motion
%     gama_diff = gama_d - gama;
%     vadali_control = rho*(eye(4,4)-(D')*(inv(D*D' + alpha*eye(3,3))*D))*gama_diff;
% u(:,i) = D'*inv(D*D')*(-Lr) + vadali_control;
%     vadali_control = rho*(eye(4,4)-(A')*(inv(A*A' + alpha*eye(3,3))*A))*gama_diff;
%     u(:,i) = A'*inv(A*A')*(-Lr(:,i)) + vadali_control;


    % u(:,i) = A*[1;0;0];
%     wheelAcc_d(:,i) = deg2rad([sin(0.02*t(i)); cos(0.02*t(i)); -cos(0.03*t(i)); -cos(0.03*t(i))]);
    
%     counter = counter + i;
%     for counter == 
%     wheelAcc_d(:,i) = u(1:4,i); 
    gamaDot_d(:,i) = u(:,i);
%         u_s(:,i) = I_w_s * (wheelAcc_d(:,i) + X(11:14)'.*omega_t);
%     end
    
   
    
%     gamaDot_d_curr = deg2rad([sind(0.02*t(i)); cosd(0.02*t(i)); -cosd(0.03*t(i)); -cosd(0.03*t(i))]);
%     gamaDot_d = deg2rad([sin(0.02*t(i)); cos(0.02*t(i)); -cos(0.03*t(i)); -cos(0.03*t(i))]);


%     gamaDot_d_next = deg2rad([sind(0.02*t(i+1)); cosd(0.02*t(i+1)); -cosd(0.03*t(i+1)); -cosd(0.03*t(i+1))]);
%     gamaDotDot_d = (gamaDot_d_next - gamaDot_d_curr)/delT;
    gamaDotDot = -K_gamaDot .* (gamaDotCurr - gamaDot_d(:,i));% + gamaDotDot_d;
    
    u_g(:,i) = Jg_g*(gamaDotDot) - ((Js_g - Jt_g)*(omega_s.*omega_t)) - (I_w_s*wheelSpeedCurr.*omega_t);
%     u_g(:,i) = u(:,i);
    
    k1 = delT*return_Steyn_CMG_Dyn(X(:,i),u_g(:,i),L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G_0,initGama,wheelSpeedCurr,I_w_t,h_dot,A);
    k2 = delT*return_Steyn_CMG_Dyn(X(:,i)+(k1/2), u_g(:,i),L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G_0,initGama,wheelSpeedCurr,I_w_t,h_dot,A);                    
    k3 = delT*return_Steyn_CMG_Dyn(X(:,i)+(k2/2), u_g(:,i),L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G_0,initGama,wheelSpeedCurr,I_w_t,h_dot,A); 
    k4 = delT*return_Steyn_CMG_Dyn(X(:,i)+k3, u_g(:,i),L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G_0,initGama,wheelSpeedCurr,I_w_t,h_dot,A);

%     [Xn_dot,H_n(:,i),T(i),T_dot(i),hs(:,i+1)] = return_X_dot_vscmg(X(:,i),u(:,i),Js_g, Jt_g, Jg_g,u_s, u_g,L_ext, N,I_sc, I_w_s,delT,G_0,initGama,I_gs,I_rw);
    [Xn_dot] = return_Steyn_CMG_Dyn(X(:,i),u_g(:,i),L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G_0,initGama,wheelSpeedCurr,I_w_t,h_dot,A);

%     if t(i) == 150
% %         H_n_30 = H_n
% %         T_30 = T(i)
% %         T_dot_40 = T_dot %(T(i) - T(i-1))/delT
%         sigma_bn_30 = X(1:3,i);
%         omega_bn_30 = X(4:6,i);
%         gama_30 = X(7:10,i);
%         wheelSpeed_30 = X(15:18,i);
%     end
    
    X_dot(:,1) = (1/6)*(k1+2*k2+2*k3+k4);
    X(:,i+1) = X(:,i) + X_dot(:,1);   
    
    if norm(X(1:3,i+1)) > 1
        X(1:3,i+1) = -(X(1:3,i+1)/(norm(X(1:3,i+1))^2));
    end

end

getSteynPlots(t,X,sigmaErr,u_g,Lr)

