function [X] = tradeoff_comrad( rou,H,Y,power,X_arbi )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here% [1]﻿Toward Dual-functional Radar-Communication Systems: Optimal Waveform Design 
[N,~] = size(H);
[~,L] = size(Y);
p1 = sqrt(rou);
p2 = sqrt(1-rou);
A = [p1*H,p2*eye(N)].'; % （17） of  [1]
B = [p1*sqrt(power)*Y;p2*X_arbi];% （17） of  [1]
Q = A'*A;% two lines above (20） of  [1]
G = A'*B;% two lines above (20） of  [1]
%G_inv = inv(G*G');
e = eig(Q);
[V,~] = eig(Q);
l_min = min(e);
for i = 1:N
    tr_fe(i) = real(trace (G*G'*V(:,i)*V(:,i)'));
end
func = @(xx)sum(tr_fe'./(e+xx*ones(N,1)))+xx*power*L; %% xx =lamada
LB = -l_min;
UB = 1;
EPSILON = 10^-4;
l = LineSearchGoldenSection(func,LB,UB,EPSILON);
X = pinv(Q+l*eye(N))*G;% （23） of  [1]
end

