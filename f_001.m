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
clear all,close all,clc;yalmip('clear');

A=diag([-1,-2]);
Ad=eye(2);

nx=2;

eps1=1e-5; 
% P=sdpvar(nx,nx,'symmetric');
% S=sdpvar(nx,nx,'symmetric');
P=sdpvar(nx,nx,'diagonal');
S=sdpvar(nx,nx,'diagonal');

F=[];
F=[F;P>=eps1*eye(nx)];
F=[F;[[P*A]+A'*P+S,[P*Ad];Ad'*P,-S]<=0*eye(2*nx)];
F=[F;0<=vec(P)<=10000];
F=[F;0<=vec(S)<=10000];

ops = sdpsettings('solver','sdpt3');
sol = optimize(F,[],ops);
sol.info

P0=value(P)
eig(P0)
S0=value(S)
eig(S0)
max(eig([[P0*A]+[P0*A]'+S0,[P0*Ad];[P0*Ad]',-S0]))
%%

%% matlab delayed-diff system simulation
clear all,close all,clc;yalmip('clear');
% BVP: bvp4c,bvp5c,bvpinit
% ODE: ode45
% DDE: dde23,ddesd,ddensd
lags=[0.2];
t_vec=[0:1e-4:20]';
fig1=figure(1); fig1.Color=[1,1,1];
for ii=1:1:5
    sol1=dde23(@dde_func,lags,@x_history,t_vec);
    x_vec = deval(sol1,t_vec);
    x1_trajectory=x_vec(1,:);
    x2_trajectory=x_vec(2,:);
    plot(x1_trajectory,x2_trajectory,...
        'LineStyle','-',...
        'LineWidth',[3],...
        'Color','r'); hold on; xlabel('x1'); ylabel('x2');
%     plot(t_vec,x1_trajectory,'LineStyle','-','LineWidth',[3],'Color','r'); hold on;
%     plot(t_vec,x2_trajectory,'LineStyle','-','LineWidth',[3],'Color','b');
%     legend('x1','x2'); xlabel('time'); ylabel('x(t)');
end

function xdot=dde_func(t,x,x_delayed)
x1=x(1);
x2=x(2);
xdot=zeros(2,1);
x1_delayed=x_delayed(1);
x2_delayed=x_delayed(2);
    A=diag([-1,-2]);
    Ad=eye(2);   
    xdot=A*[x1;x2]+Ad*[x1_delayed;x2_delayed];
end

function x=x_history(t)
%     x=ones(2,1);
    x=0.1*ones(2,1)+rand(2,1)*0.3;
%     x=[sin(t);sin(2*t)];
end









%