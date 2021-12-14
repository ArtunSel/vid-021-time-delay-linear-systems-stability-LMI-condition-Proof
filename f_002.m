clear all,close all,clc;
%%
% this is for SDPT3
current_dir=pwd;
cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);

% this is for YALMIP
addpath(genpath('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\YALMIP-master'))
yalmip('clear');
%% BEFORE RUNNING THE CODE ADD "YALMIP" AND "SDPT3" libraries
%% yalmip analysis
%     set(findall(gcf,'type','line'),'linewidth',[1]);
clear all,close all,clc;yalmip('clear');
A=[-2,0,1;0,-3,0;1,0,-2];
Ad=[-1,1,1;2,-1,1;0,0,-1];

nx=3;           % state-dimension
dmax=0.3;       % max delay [max value 0.3 seconds]
eps1=1e-3;      % epsilon for numerical accuracy
% define the decision variables
X=sdpvar(nx,nx,'symmetric');
beta=sdpvar(1,1,'symmetric');
% enter the constraints
F=[];
F=[F;X>=eps1*eye(nx)];
M_11=[(A+Ad)*X]+[(A+Ad)*X]'+dmax*Ad*Ad';
M_12=dmax*X*A';
M_13=dmax*X*Ad';
M_21=dmax*A*X;
M_22=-dmax*beta*eye(nx);
M_23=zeros(nx);
M_31=dmax*Ad*X;
M_32=zeros(nx);
M_33=-dmax*(1-beta)*eye(nx);
M=[M_11,M_12,M_13;M_21,M_22,M_23;M_31,M_32,M_33];
F=[F;M<=-eps1*eye(3*nx)];
F=[F;-1e2<=vec(X)<=1e2];
F=[F;eps1<=beta<=1-eps1];
% solve the opt-problem
ops = sdpsettings('solver','sdpt3');
sol = optimize(F,[],ops);
sol.info
% get the decision variables
X=value(X)
P=inv(X)
beta=value(beta)
%% check if there is any error
X=value(X)
P=inv(X)
beta=value(beta)

eig(P)
M_11=[(A+Ad)*X]+[(A+Ad)*X]'+dmax*Ad*Ad';
M_12=dmax*X*A';
M_13=dmax*X*Ad';

M_21=dmax*A*X;
M_22=-dmax*beta*eye(nx);
M_23=zeros(nx);

M_31=dmax*Ad*X;
M_32=zeros(nx);
M_33=-dmax*(1-beta)*eye(nx);
M=[M_11,M_12,M_13;M_21,M_22,M_23;M_31,M_32,M_33];
eig(P)
eig(M)
%%

%% matlab delayed-diff system simulation
clear all,close all,clc;yalmip('clear');

lags=[0.2];
t_vec=[0:1e-4:20]';
fig1=figure(1);
fig1.Color=[1,1,1];
for ii=1:1:5
    sol1=dde23(@dde_func,lags,@x_history,t_vec);
    x_vec = deval(sol1,t_vec);
    x1_trajectory=x_vec(1,:);
    x2_trajectory=x_vec(2,:);
    x3_trajectory=x_vec(3,:);
    plot3(x1_trajectory,x2_trajectory,x3_trajectory,...
        'LineStyle','-',...
        'LineWidth',[1],...
        'Color','r'); hold on;    
%     plot(t_vec,x1_trajectory,'LineStyle','-','LineWidth',[1],'Color','r'); hold on;
%     plot(t_vec,x2_trajectory,'LineStyle','-','LineWidth',[1],'Color','g');
%     plot(t_vec,x3_trajectory,'LineStyle','-','LineWidth',[1],'Color','b');
%     legend('x1','x2','x3');xlabel('time');ylabel('x(t)');
end
function xdot=dde_func(t,x,x_delayed)
x1=x(1);
x2=x(2);
x3=x(3);
xdot=zeros(3,1);
x1_delayed=x_delayed(1);
x2_delayed=x_delayed(2);
x3_delayed=x_delayed(3);
A=[-2,0,1;0,-3,0;1,0,-2];
Ad=[-1,1,1;2,-1,1;0,0,-1];
    xdot=A*[x1;x2;x3]+Ad*[x1_delayed;x2_delayed;x3_delayed];
end
function x=x_history(t)
%     x=ones(2,1);
    x=0.1*ones(3,1)+rand(3,1)*0.3;
end









%