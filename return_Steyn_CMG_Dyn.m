function [Xn_dot] = return_Steyn_CMG_Dyn(X, u, L_ext, N,I_sc, I_w_s,Js_g, Jt_g, Jg_g,h,G,initGama,wheelSpeedCurr,I_w_t,h_dot,A)

sigma_bn = X(1:3);
omega_bn = X(4:6);
omega_tilde_bn = [0, -omega_bn(3), omega_bn(2); ...
                  omega_bn(3), 0, -omega_bn(1); ...
                  -omega_bn(2), omega_bn(1), 0];
gamaCurr = X(7:10);
gamaDotCurr = X(11:14);
% wheelSpeedCurr = X(15:18);
              
J_g = diag([I_w_s,I_w_t,I_w_t]);   
% 
gs_0 = G(1:3,1:4);
gt_0 = G(4:6,1:4);
gg   = G(7:9,1:4);
gama = gamaCurr;
gamaDot = gamaDotCurr;

for k = 1:N
    delAngle(k) = gama(k) - initGama(k); % subtracting from initial gama instead of previous gama because error builds up using previous gama instead of the absolute or constant value of initial gama. also this euation wouldnt work for previous gama. 
    gs(:,k) = [cos(delAngle(k))*gs_0(1:3,k) + sin(delAngle(k))*gt_0(1:3,k)];
    gt(:,k) = [-sin(delAngle(k)).*gs_0(1:3,k) + cos(delAngle(k)).*gt_0(1:3,k)];
end
    
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


I = I_sc + J_b_1 + J_b_2 + J_b_3 + J_b_4;

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


f_sigma = (BmatMRP(sigma_bn)*omega_bn)/4;
f_gama = gamaDotCurr; % equals to gamaDot
% f_gamaDot(k,:) = u(k
total = 0;

for k = 1:N
%     f_omega_wheels(k,:) = us(k) - (I_w_s*f_gama(k)*omega_t(k)); % equals to omegaDot
%     f_gamaDot(k,:) = u(k) + ((Js_g - Jt_g)*omega_s(k)*omega_t(k)) + I_w_s*wheelSpeedCurr(k)*omega_t(k);
    f_gamaDot(k,:) = u(k) + ((I_w_s - I_w_t)*omega_s(k)*omega_t(k)) + I_w_s*wheelSpeedCurr(k)*omega_t(k);

%     total = total + (gs(1:3,k)*(Js_g*f_gama(k)*omega_t(k) - (Jt_g-Jg_g)*omega_t(k)*f_gama(k)) + gt(1:3,k)*((Js_g*omega_s(k) + I_w_s*wheelSpeedCurr(k))*f_gama(k) - (Jt_g + Jg_g)*omega_s(k)*f_gama(k) + I_w_s*wheelSpeedCurr(k)*omega_g(k)) - gg(1:3,k)*I_w_s*wheelSpeedCurr(k)*omega_t(k));
%     total = total + (gs(1:3,k)*(I_w_s*f_gama(k)*omega_t(k)) + gt(1:3,k)*(I_w_s*f_gama(k)*(omega_s(k) + wheelSpeedCurr(k)) - (2*I_w_t)*omega_s(k)*f_gama(k) + I_w_s*wheelSpeedCurr(k)*omega_g(k)) - gg(1:3,k)*I_w_s*wheelSpeedCurr(k)*omega_t(k));
    total = total + (gs(1:3,k)*I_w_s*f_gama(k)*omega_t(k));
end

% hs = I_w_s*[(omega_s_1+wheelSpeedCurr(1));(omega_s_2+wheelSpeedCurr(2));(omega_s_3+wheelSpeedCurr(3));(omega_s_4+wheelSpeedCurr(4))];

% f_omega_sc = -(omega_tilde_bn*I*omega_bn) - total;
f_omega_sc = -(omega_tilde_bn*I*omega_bn) - A*gamaDotCurr - omega_tilde_bn*h;
% f_omega_sc = -(omega_tilde_bn*I*omega_bn) - omega_tilde_bn*h;

M = [  eye(3),        zeros(3,3),    zeros(3,N),    zeros(3,N); ...
        zeros(3,3),    I,         zeros(3,N),    I_w_t*gg; ...
        zeros(N,3),    zeros(N,3),   eye(N,N),      zeros(N,N); ...
        zeros(N,3),    I_w_t*gg',    zeros(N,N),    I_w_t*eye(N,N)];    

Xn_dot = inv(M)*[f_sigma; f_omega_sc; f_gama; f_gamaDot];

% BN = MRP2C(sigma_bn);
% NB = BN';

% omega_gb = [omega_s,omega_t,omega_g]; % omega_s_1_g and so on
% omega_nbg = NB*BG*omega_g; % omega_s_1_b and so on
% H_bn = NB*I_sc*BN*(NB*omega_bn);
% H_bn = I_sc*(omega_bn);
% 
% 
% total_e_wg = 0;
% total_t_dot = 0;
% for m = 1:N
%     H_wn_gn(:,m) = [(Js_g*omega_gb(m,1) + I_w_s*wheelSpeedCurr(m));(Jt_g*omega_gb(m,2));(Jg_g*omega_gb(m,3) + Jg_g*gamaDot(m))];
%     total_e_wg = total_e_wg + 0.5*(I_w_s*(wheelSpeedCurr(m)+omega_s(m))^2 + (I_gs*omega_s(m)^2) + Jt_g*omega_t(m)^2 + Jg_g*(omega_g(m) + gamaDot(m))^2);
%     total_t_dot = total_t_dot + gamaDot(m)*u_g(m) + wheelSpeedCurr(m)*u_s(m);
% end
% 
% T = (0.5*omega_bn'*I_sc*omega_bn) + total_e_wg;
% T_dot = omega_bn'*L_ext + total_t_dot;
% 
% H_n = (NB*H_bn)+(NB*BG_1*H_wn_gn(:,1)) +  (NB*BG_2*H_wn_gn(:,2)) + (NB*BG_3*H_wn_gn(:,3)) + (NB*BG_4*H_wn_gn(:,4));
% 
end