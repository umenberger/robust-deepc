function res = cdc_plot(ops_deepc,varargin)

sys_1 = varargin{1};
sys_2 = varargin{2};
sys_3 = varargin{3};

Tini = ops_deepc.Tini;

Tcs = sys_1.t - 1 - Tini;

[~,u1,w1,y1] = sys_1.get_history();

w1 = w1(:,Tini+1:end);

u1 = u1(:,Tini+1:end);
y1 = y1(:,Tini+1:end);

[~,u2,w2,y2] = sys_2.get_history();

w2 = w2(:,Tini+1:end);

u2 = u2(:,Tini+1:end);
y2 = y2(:,Tini+1:end);

[~,u3,w3,y3] = sys_3.get_history();

w3 = w3(:,Tini+1:end);

u3 = u3(:,Tini+1:end);
y3 = y3(:,Tini+1:end);

r = ops_deepc.reference;

umin = ops_deepc.umin;
umax = ops_deepc.umax;

if isprop(ops_deepc,'ymax')
   
    ymax = ops_deepc.ymax;
    ymin = ops_deepc.ymin;
    
    output_constraints = 1;
    
else
    
    output_constraints = 0;
    
end

resolve_period = ops_deepc.resolve_period;

nu = size(u2,1);
ny = size(y2,1);

lw = 2;

plot_w = 0;

title_on = 0;

if plot_w
    nsp = 3;
else
    nsp = 2;
end

fs = 15;

L = 3;

% colormap from linspecer
if exist('linspecer','file')
    colors = [linspecer(L), 0.7*ones(L,1)];
else
    tmp = [    0.3467    0.5360    0.6907,
               0.9153    0.2816    0.2878,
               0.4416    0.7490    0.4322];

    colors = [tmp 0.7*ones(L,1)];
end


figure

ax = subplot(nsp,1,1);

for i = 1:ny  
    
    rmax = max(r(i,1:Tcs));
    rmin = min(r(i,1:Tcs));
    
    plot(1:Tcs,r(i,1:Tcs),'linewidth',4,'color',[0.5,0.5,0.5,0.5])
    hold on
    plot(1:Tcs,y1(i,1:Tcs),'color',colors(1,:),'linewidth',lw)
    plot(1:Tcs,y2(i,1:Tcs),'color',colors(2,:),'linewidth',lw) 
    plot(1:Tcs,y3(i,1:Tcs),'color',colors(3,:),'linewidth',lw)
    
    if output_constraints
        
        plot(1:Tcs,ymax(i,1:Tcs),'k-.','linewidth',lw)
%         plot(1:Tcs,ymin(i,1:Tcs),'k')
        
    end
    
    
    for j = 1:ceil(Tcs/resolve_period)
       plot((1+(j-1)*resolve_period)*[1,1],[min(ops_deepc.ymin(i,:)),max(ops_deepc.ymax(i,:))],'linewidth',0.8,'color',[0,0,0,0.5])  
    end   
    
    
    set(gca,'ylim',[min(r(i,:))-0.2,1.1*max(ops_deepc.ymax(i,:))])
     
end


ylabel('output','interpreter','latex','fontsize',fs)
if title_on
    title(['g reg = ' num2str(ops_deepc.soln_regularization) ', '...
           'OL-inf = ' num2str(tracking_error(y1,r,inf)) ', CL-inf = ' num2str(tracking_error(y2,r,inf)) ...
           ', OL-1 = ' num2str(tracking_error(y1,r,1)) ', CL-1 = ' num2str(tracking_error(y2,r,1))])
end

methods = {'reference','open-loop','closed-loop','non-robust','constraint','re-compute'};

legend(methods,'interpreter','latex','fontsize',fs)
legend boxoff

p = get(ax,'pos');
p(2) = 0.8*p(2);
p(4) = 1.4*p(4);
set(ax,'pos',p)

% set(gca,'fontsize',fs)

ax = subplot(nsp,1,2);

plot(nan,nan,'b','linewidth',lw)
hold on
plot(nan,nan,'r--','linewidth',lw)
plot(1:Tcs,umax*ones(1,Tcs),'k-.','linewidth',1)
plot(1:Tcs,umin*ones(1,Tcs),'k-.','linewidth',1)
for i = 1:nu  
    plot(1:Tcs,u1(i,1:Tcs),'color',colors(1,:),'linewidth',lw)
    hold on
    plot(1:Tcs,u2(i,1:Tcs),'color',colors(2,:),'linewidth',lw)   
    plot(1:Tcs,u3(i,1:Tcs),'color',colors(3,:),'linewidth',lw)
    
    for j = 1:ceil(Tcs/resolve_period)
       plot((1+(j-1)*resolve_period)*[1,1],[umin,umax],'linewidth',0.8,'color',[0,0,0,0.5]) 
    end      
    
end
xlabel('time','interpreter','latex','fontsize',fs)
ylabel('input','interpreter','latex','fontsize',fs)
set(gca,'ylim',1.1*[umin,umax])


p = get(ax,'pos');
p(4) = 0.8*p(4);
set(ax,'pos',p)

if plot_w
    subplot(nsp,1,3)

    plot(1:Tcs,w1)
    hold on
    plot(1:Tcs,w2,'r--')
    xlabel('time')
    ylabel('disturbance')
    legend('disturbance')
end


% errors
res.norm_1 = [tracking_error(y1,r,1); tracking_error(y2,r,1); tracking_error(y3,r,1)];

% check for constraint violation

if sum(y3 > ops_deepc.ymax(1:length(y3))) || sum(y3 < ops_deepc.ymin(1:length(y3)))
    violation = 1;
else
    violation = 0;
end

res.violation = violation;




end