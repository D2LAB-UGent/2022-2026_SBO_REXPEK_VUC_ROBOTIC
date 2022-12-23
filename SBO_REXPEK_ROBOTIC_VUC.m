function SBO_REXPEK_ROBOTIC_VUC

clc
close all
clear all

test = false;

answer = {};

% get some input
if ~test
    while isempty(answer)
        prompt = {'Your name:'};
        dlgtitle = 'User information';
        dims = [1 35];
        definput = {'',''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer)
            answer = {'anonymous_user'};
        end
        subjectName = answer{1};
    end
end

% do the actual survey
seq = do_survey; close all;

% send email
if ~test
    answer = questdlg('Would you like to share your data?', ...
        'Data processing', ...
        'Yes.','No.','Yes.');
    switch answer
        case 'Yes.'
            tok = clock;
            token = '';
            for j = 2:length(tok)-1
                token = [token,num2str(tok(j))];
            end
            token = [token,num2str(ceil(tok(6)))];
            fileName = [subjectName,'_',token,'.csv'];
            fileID = fopen(fileName,'at');
            fprintf(fileID,['i',', ','a',', ','A',', ','N',', ','T',', ','D',', ','F','\n']);
            for i = 1:size(seq,1)
                fprintf(fileID,[num2str(i),', ',...
                    num2str(seq(i,1)),', ',...
                    num2str(seq(i,2)),', ',...
                    num2str(seq(i,3)),', ',...
                    num2str(seq(i,4)),'\n']);
            end
            fclose(fileID);
            setpref('Internet','SMTP_Server','smtprelay.UGent.be');
            setpref('Internet','E_mail','tom.lefebvre@ugent.be');
            sendmail('tom.lefebvre@ugent.be',...
                '2022-2026 SBO REXPEK : VIRTUAL USE CASE : ROBOT : DATA',...
                [subjectName,' send this.'],...
                [pwd,'\',fileName]);
            delete([pwd,'\',fileName]);
        case 'No.'
    end
end

clear all;

end

% gui related functions

function do = safe_prompt(main_prompt)
do = false;
ver = input(['Do you want to ',main_prompt,'? [y/n] '],'s');
if ver == 'y'
    do = true;
end
disp(newline);
end

function out = data_prompt(main_prompt,verification_prompt)
do = true;
while do
    out = input(main_prompt,'s');
    ver = input([verification_prompt,out,'? [y/n] '],'s');
    if ver == 'y'
        do = false;
    end
end
disp(newline);
end

function output = do_survey

% gui

screen = get(0,'screensize');
main = uifigure('Name','SBO REXPEK 2022-2026 : Robotic virtual use-case. UGent-D2LAB/FM-MIRO ©.','Position',[screen(3)/2-690 screen(4)/2-398 1380 796]);

pnl1 = uipanel(main,'Position',[20 20 340 756],'title','Main','BackgroundColor','w');
pnl2 = uipanel(main,'Position',[380 20 980 718],'title','Experiment','BackgroundColor','w');

pnl1_1 = uipanel(pnl1,'Position',[20 498 300 219],'title','Control panel');
pnl1_2 = uipanel(pnl1,'Position',[20 259 300 219],'title','Parameter setting');
pnl1_3 = uipanel(pnl1,'Position',[20 20 300 219],'title','History');

sld_1 = uislider(pnl1_2,'Position',[50 149 200 3],'limits',[0 1],'value',0);
sld_2 = uislider(pnl1_2,'Position',[50 98.25 200 3],'limits',[0 2],'value',1);
sld_3 = uislider(pnl1_2,'Position',[50 47.75 200 3],'limits',[2 20],'value',15);

uibutton(pnl1_1,'text','run experiment','position',[100 40 100 40],'ButtonPushedFcn',@(src,event) run_experiment);
uibutton(pnl1_1,'text','stop tuning','position',[100 120 100 40],'ButtonPushedFcn',@(src,event) close(main));

pnl2_1 = uipanel(pnl2,'Position',[20 20 300 658],'title','Visualization');
pnl2_2 = uipanel(pnl2,'Position',[340 359 300 319],'title','Inputs');
pnl2_3 = uipanel(pnl2,'Position',[660 359 300 319],'title','Outputs');
pnl2_4 = uipanel(pnl2,'Position',[340 20 620 319],'title','Animation');

uibutton(pnl2_4,'text','run animation','position',[140 239 100 40],'ButtonPushedFcn',@(src,event) run_animation);
uibutton(pnl2_4,'text','stop animation','position',[380 239 100 40],'ButtonPushedFcn',@(src,event) stop_animation);

% create global variables

global animate_flag ax2_1 ax2_2 ax2_3 ax2_4
animate_flag = false;

output = zeros(0,4);

It = [];
Xt = [];
Ut = [];
Ot = [];
context = [];

waitfor(main);

    function run_experiment
        
        if ~isempty(ax2_1) %#ok<*NODEF>
            cla(ax2_1); drawnow
            cla(ax2_2); drawnow
            cla(ax2_3); drawnow
        end
        
        if ~isempty(ax2_4) %#ok<*NODEF>
            cla(ax2_4); drawnow
        end
        
        animate_flag = false;
        
        val1 = get(sld_1,'value');
        val2 = get(sld_2,'value');
        val3 = get(sld_3,'value');
        
        [result,context,It,Xt,Ut,Ot] = experiment(val1,val2,val3);
        
        output = [output; val1 val2 val3 result];
        
        ax1_1 = uiaxes(pnl1_3,'Position',[10 10 135 180],'box','on');
        ax1_2 = uiaxes(pnl1_3,'Position',[155 10 135 180],'box','on');
        
        plot_history(ax1_1,ax1_2,output);
        
        ax2_1 = uiaxes(pnl2_1,'Position',[10 10 280 619],'box','on'); hold(ax2_1,'on');
        ax2_2 = uiaxes(pnl2_2,'Position',[10 10 280 280],'box','on'); hold(ax2_2,'on');
        ax2_3 = uiaxes(pnl2_3,'Position',[10 10 280 280],'box','on'); hold(ax2_3,'on');
        
        plot_results(ax2_1,ax2_2,ax2_3,context,It,Xt,Ut,Ot);
        
    end

    function run_animation
        if ~isempty(Xt)
            window = [-10 10 -5 25 0 1];
            ax2_4 = uiaxes(pnl2_4,'Position',[10 10 600 210],'box','on','view',[-75 30]);
            camproj(ax2_4,'perspective'); daspect(ax2_4,[1 1 1]); axis(ax2_4,window);
            animate_flag = true;
            aniTrajOt(ax2_4,context,Xt,Ot,[.75 .75 .75],window);
        end
    end

    function stop_animation
        animate_flag = false;
    end

end

function plot_history(ax1,ax2,output)
h1 = plot(ax1,1:size(output,1),output(:,1:3),'linewidth',2);
legend(h1,{'$\alpha$','$A_0$','$N$'},'interpreter','latex');
h2 = plot(ax2,1:size(output,1),output(:,4),'linewidth',2);
legend(h2,{'$\int {d}x$'},'interpreter','latex'); set(ax2,'YScale','log')
end

function plot_results(ax1,ax2,ax3,context,It,Xt,Ut,Ot)

for i = 1:length(It)-1
    switch Ot(It(i))
        case 1
            visTraj(ax1,Xt(:,It(i):It(i+1)-1),context.geom1,'g');
        case 2
            visTraj(ax1,Xt(:,It(i):It(i+1)-1),context.geom2,'r');
    end
end

h1 = plot(ax2,[Xt;sqrt((Xt(1,:)-context.xg).^2+(Xt(2,:)-context.yg).^2)]','linewidth',1); xlim(ax2,[0 length(Xt)-1]);
plot(ax2,[0 length(Xt)],[context.xg context.xg],'k--','linewidth',1);
plot(ax2,[0 length(Xt)],[context.yg context.yg],'k--','linewidth',1);
legend(h1,{'$x$','$y$','$\theta$','$c$','$|d|$'},'interpreter','latex','location','northeast');

h2 = plot(ax3,Ut','linewidth',1); xlim(ax3,[0 length(Xt)-1]);
legend(h2,{'$u$','$v$'},'interpreter','latex','location','northeast');

end

% experiment

function [val,context,It,Xt,Ut,Ot] = experiment(a,b,c)

% define context

context.type = 1;

context.a = 4;
context.b = 1;
context.r = .1;

context.x0 = 5;
context.y0 = 0;
context.t0 = 0;

context.xg = 0;
context.yg = 20;

% define control

control.alpha = a;
control.A0 = b;
control.N = ceil(c);

% do experiment

[context,It,Xt,Ut,Ot] = simPush(control,context);

val = sum(sqrt(sum((Xt(1:2,1:end-1)-Xt(1:2,2:end)).^2,1)));

end

% simulation functions

function [context,It,Xt,Ut,Ot] = simPush(control,context)

type = context.type;

a = context.a;
b = context.b;
r = context.r;

geom1.a = a;
geom1.b = b;
geom1.r = r;

[~,~,~,~,f1] = dynamics(geom1,type);

geom2.a = b;
geom2.b = a;
geom2.r = r;

context.geom1 = geom1;
context.geom2 = geom2;

[~,~,~,~,f2] = dynamics(geom2,type);

x = context.x0;
y = context.y0;
t = context.t0;

xg = context.xg;
yg = context.yg;

alpha = control.alpha;
A0 = control.A0;
N = control.N;

orientation = 1;

[tn,cn,orientation] = update(x,y,t,a,b,xg,yg,orientation);

It = 1;
Xt = [];
Ut = [];
Ot = [];

k = 0;

while sqrt((x-xg)^2+(y-yg)^2) > 1 && k <= 1e3
    
    k = k+1;
    
    utref = [   ones(1,N);
        zeros(1,N)];
    theta = (1-alpha)*pi2pi(atan2(yg-y,xg-x))+alpha*pi2pi(tn+pi/2);
    utref = rot(theta)*utref;
    xt = zeros(4,N);
    xt(:,1) = [x;y;tn;cn];
    ut = zeros(2,N);
    ot = orientation*ones(1,N);
    
    A = A0*10e-3*sqrt((x-xg)^2+(y-yg)^2);
    
    switch orientation
        case 1
            for t = 1:N-1
                ut(:,t) = A*rot(xt(3,t))'*utref(:,t);
                if abs(xt(4,t)) > geom1.a/2
                    xt = xt(:,1:t);
                    ut = ut(:,1:t);
                    ot = ot(:,1:t);
                    break
                end
                xt(:,t+1) = xt(:,t) + f1(xt(:,t),ut(:,t));
                if t == N-1
                    ut(:,t+1) = ut(:,t);
                end
            end
        case 2
            for t = 1:N-1
                ut(:,t) = A*rot(xt(3,t))'*utref(:,t);
                if abs(xt(4,t)) > geom2.a/2
                    xt = xt(:,1:t);
                    ut = ut(:,1:t);
                    ot = ot(:,1:t);
                    break
                end
                xt(:,t+1) = xt(:,t) + f2(xt(:,t),ut(:,t));
                if t == N-1
                    ut(:,t+1) = ut(:,t);
                end
            end
    end
    
    It = [It It(end)+size(xt,2)];
    Xt = [Xt xt];
    Ut = [Ut ut];
    Ot = [Ot ot];
    
    x = xt(1,end);
    y = xt(2,end);
    t = xt(3,end);
    
    switch orientation
        case 1
            [tn,cn,temp] = update(x,y,t,geom1.a,geom1.b,xg,yg,orientation);
        case 2
            [tn,cn,temp] = update(x,y,t,geom2.a,geom2.b,xg,yg,orientation);
    end
    orientation = temp;
    
end

end

function [t,c,orientation] = update(x,y,t,a,b,xg,yg,orientation)
theta = atan2(yg-y,xg-x);
alpha = atan2(b,a);
beta = pi/2-alpha;
dt = pi2pi(theta-t-pi/2);
if dt < -beta-2*alpha       % AREA I
    change = 0;
    t = pi2pi(t + pi);
    c = b/2*tan(dt+pi);
elseif dt < -beta           % AREA II
    change = 1;
    t = pi2pi(t-pi/2);
    c = a/2*tan(pi/2+dt);
elseif dt < beta            % AREA III
    change = 0;
    c = b/2*tan(dt);
elseif dt < beta+2*alpha    % AREA IV
    change = 1;
    t = pi2pi(t+pi/2);
    c = -a/2*tan(pi/2-dt);
elseif dt <= pi             % AREA V
    change = 0;
    t = pi2pi(t + pi);
    c = b/2*tan(dt);
end
if change
    switch orientation
        case 1
            temp = 2;
        case 2
            temp = 1;
    end
    orientation = temp;
end
end

function t = pi2pi(t)
if t < -pi
    t = t+2*pi;
elseif t > pi
    t = t-2*pi;
end
end

function [dx,dy,dt,dc,f] = dynamics(geom,type)
a = geom.a;
b = geom.b;
r = geom.r;
if nargin < 2
    type = 1;
end
switch type
    case 1
        beta = (log(sqrt(a^2 + b^2) + b)*a^3 - log(a)*a^3 + log(b)*b^3 - log(sqrt(a^2 + b^2) - a)*b^3 + 2*sqrt(a^2 + b^2)*a*b)/(12*a*b);
    case 2
        beta = sqrt((a^2+b^2)/12);
end
dx = @(c,u,v) 0;
dy = @(c,u,v) beta^2/(beta^2 + c^2)*v;
dt = @(c,u,v) c*v/(beta^2 + c^2);
dc = @(c,u,v) u - (b/2+r)*c/(beta^2 + c^2)*v;

f = @(x,u) [rot(x(3))*[dx(x(4),u(1),u(2));
    dy(x(4),u(1),u(2))];
    dt(x(4),u(1),u(2));
    dc(x(4),u(1),u(2))];
end

% visualization functions

function aniTrajOt(ax,context,Xt,Ot,color,window)

global animate_flag;

r = context.geom1.r;

xg = cos(linspace(0,2*pi));
yg = sin(linspace(0,2*pi));

[xc,yc,zc] = cylinder; xc = r*xc; yc = r*yc; zc = .5*zc;

T = size(Xt,2);
dd = zeros(2,T);
for t = 1:T
    switch Ot(t)
        case 1
            dd(:,t) = rot(Xt(3,t))*[Xt(4,t);-context.geom1.b/2-context.geom1.r] + Xt(1:2,t);
        case 2
            dd(:,t) = rot(Xt(3,t))*[Xt(4,t);-context.geom2.b/2-context.geom2.r] + Xt(1:2,t);
    end
end

for t = 1:size(Xt,2)
    if animate_flag
        cla(ax);
        daspect(ax,[1 1 1]); camproj(ax,'perspective'); view(ax,[-75 30]); axis(ax,window);
        fill(ax,context.xg+xg,context.yg+yg,'r')
        switch Ot(t)
            case 1
                vis3(ax,Xt(1,t),Xt(2,t),Xt(3,t),context.geom1.a,context.geom1.b,color); hold(ax,'on');
                surf(ax,xc+dd(1,t),yc+dd(2,t),zc,'edgecolor','none');
            case 2
                vis3(ax,Xt(1,t),Xt(2,t),Xt(3,t),context.geom2.a,context.geom2.b,color); hold(ax,'on');
                surf(ax,xc+dd(1,t),yc+dd(2,t),zc,'edgecolor','none');
        end
        drawnow
    else
        break
    end
end

pause(1);

end

function vis3(ax,x,y,theta,a,b,color)
height = .25;
R = rot(theta);
p = diag([a/2,b/2])*[1 1 -1 -1 1;-1 1 1 -1 -1];
pp = R*p + [x;y];
fill3(ax,pp(1,:),pp(2,:),zeros(size(pp(1,:))),color,'edgecolor','k');
fill3(ax,pp(1,:),pp(2,:),0.25*ones(size(pp(1,:))),color,'edgecolor','k');
z = height*[1 1 0 0 1];
for i = 1:4
    fill3(ax,pp(1,[0 1 1 0 0]+i),pp(2,[0 1 1 0 0]+i),z,color,'edgecolor','k');
end
end

function R = rot(theta)
R = [[cos(theta) -sin(theta)];[sin(theta) cos(theta)]];
end

function visTraj(ax,xt,geom,color)
a = geom.a;
b = geom.b;
r = geom.r;
if nargin < 3
    color = 'g';
end
circ.x = r*cos(linspace(0,2*pi,20));
circ.y = r*sin(linspace(0,2*pi,20));
[~,T] = size(xt);
dd = zeros(2,T);
for t = 1:T
    vis(ax,xt(1,t),xt(2,t),xt(3,t),a,b,color);
    dd(:,t) = rot(xt(3,t))*[xt(4,t);-b/2-r] + xt(1:2,t);
    h = fill(ax,circ.x+dd(1,t),circ.y+dd(2,t),'y'); set(h,'facealpha',.5);
end
end

function vis(ax,x,y,theta,a,b,color)
R = rot(theta);
p = diag([a/2,b/2])*[1 1 -1 -1 1;-1 1 1 -1 -1];
pp = R*p + [x;y];
h = fill(ax,pp(1,:),pp(2,:),color,'edgecolor','k'); set(h,'facealpha',.1);
end