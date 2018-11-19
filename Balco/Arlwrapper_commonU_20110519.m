% Wrapper script for three-parameter least squares solution

clear all; close all;

%% data entry
% From table supplied by Paul Jan 23rd
% Added 1 yr uncertainty for Vesuvius

db = [1.9190E-03	1.00E-06	6.7890E-05	5.5000E-07
        17.42000	0.03000	0.60556	0.00249
        125.65000	0.17000	4.54940	0.00610
        133.06000	0.26000	4.82050	0.01530
        201.27000	0.13000	7.46360	0.00720
        242.14000	0.45000	9.06190	0.01220
        252.17000	0.36000	9.49830	0.01020
        252.40000	0.33000	9.49180	0.00410
        252.60000	0.33000	9.51020	0.00590
        253.02000	0.50000	9.49580	0.01220
        253.24000	0.33000	9.52950	0.00680
        253.68000	0.30000	9.53910	0.00850
        257.30000	0.30000	9.69870	0.00530
        327.99000	0.45000	12.60390	0.01480
        454.59000	0.56000	18.09000	0.05230
        1094.20000	1.15000	52.90110	0.11620
        2067.50000	1.40000	135.63890	0.18210];

data.ri = db(:,3);
data.delri = db(:,4);

data.ti = 1e6*db(:,1);
data.delti = 1e6*db(:,2);

% Other constraints

% Order: K le lb ltot
data.a = [1.6408E-03
5.8000E-11
4.8840E-10
5.5545E-10];

data.dela = [4.7000E-06
7.0000E-13
4.9000E-12
1.0900E-12];

% Do residence time correction?
doResT = 1;

if doResT == 1;
    data.ti(2:end) = data.ti(2:end) - 90000;
    data.delti(2:end) = sqrt( (data.delti(2:end)).^2 + 77000.^2);
end;

%% Try solving 

opts = optimset('fminsearch');
%opts = optimset(opts,'display','notify');
opts = optimset(opts,'MaxIter',20000);
opts = optimset(opts,'MaxFunEvals',20000);
opts = optimset(opts,'TolX',1e-6);
opts = optimset(opts,'TolFun',1e-6);

% Set up x0 
% x0 is randomly generated
% Generate with t's as input
x0 = [((randn(3,1).*data.dela(1:3))+data.a(1:3))' (randn(size(data.ti')).*0.05+1).*data.ti'];
% Generate with r's as input - must change Arls.m code as well
% x0 = [((randn(3,1).*data.dela(1:3))+data.a(1:3))' (randn(size(data.ri')).*0.05+1).*data.ri'];

% Solve
[xopt,fval,exitflag,output] = fminsearch(@(x) Arls_nol_20110519(x,data),x0,opts);

% get residuals
misses = Arls_nol_20110519(xopt,data,1);
% Report
disp(['K = ' sprintf('%0.4f',xopt(1).*1e3)]);
disp(['lambda_e = ' sprintf('%0.4f',xopt(2).*1e10)]);
disp(['lambda_b = ' sprintf('%0.4f',xopt(3).*1e10)]);
disp(['lambda_tot = ' sprintf('%0.4f',(xopt(3)+xopt(2)).*1e10)]);
disp(['fval = ' num2str(fval)]);


% Figure corresponding FCT age
tFCT = (1./(xopt(2)+xopt(3))).*log(xopt(1).*((xopt(2)+xopt(3))./xopt(2)) + 1);

disp(['tFCT = ' num2str(tFCT./1e6)]);

% Report better age estimates for samples

for a = 2:length(data.ti);
    disp(sprintf('%0.4f',xopt(3+a)./1e6));
end;

% Report better R estimates for samples

bettert = xopt(4:end);
betterR = (xopt(2)./(xopt(2)+xopt(3))).*(exp(bettert.*(xopt(2)+xopt(3)))-1)./xopt(1);

for a = 1:length(data.ti);
    disp(sprintf('%0.4f',betterR(a)));
end;

%% MC error propagation


% Must now transform t's into rhos (206*/238 ratios) for MCS.
% This allows incorporating common error in l238. 

% Define U decay constant - Jaffey value
l238 = 1.55125e-10;
dell238 = 0.00083e-10;

data.rhoi = exp(data.ti.*l238)-1;
% Also unwrap error on rhos from error on t's
dtdl = (-(log(1+data.rhoi))./(l238.^2));
dtdrho = 1./(l238.*(1+data.rhoi));
K1 = (data.delti.^2) - (dtdl.*(dell238)).^2;
data.delrhoi = sqrt(K1)./dtdrho;

% set saveFlag to 1 and set name of file to save results
saveFlag = 0;
save_fname = 'results_20110519_commonU_noLtot_4000pts';


if exist([save_fname '.mat']);
    error('Filename exists. Stopping.');
end;


numits = 4000;
results = zeros(numits,length(x0));
fval = zeros(numits,1);
opts = optimset(opts,'display','notify');

% Sample l238
l238s = randn(1,numits).*dell238 + l238;

for a = 1:numits;
    data1 = data;
    % Sample rhos
    data1.rhoi = randn(size(data.rhoi)).*data.delrhoi + data.rhoi;
    % Transform to t's
    data1.ti = (1./l238s(a)).*log(1 + data.rhoi);
    % Vesuvius is done by sampling t directly
    data1.ti(1) = randn(1,1).*data.delti(1) + data.ti(1);
    % Now sample ri's
    data1.ri = randn(size(data1.ri)).*data.delri + data.ri;
    % Sample a's
    data1.a = randn(size(data1.a)).*data.dela + data.a;
    % x0 is randomly generated
    x0 = [((randn(3,1).*data.dela(1:3))+data.a(1:3))' (randn(size(data.ti')).*0.05+1).*data.ti'];
    [thisx,thisfval] = fminsearch(@(x) Arls_nol_20110519(x,data1),x0,opts);   
    results(a,:) = thisx;
    fvals(a) = thisfval;
    if round(a./10) == a./10;
        disp(int2str(a));
    end;   
end;
    


if saveFlag;
    eval(['save ' save_fname ' results xopt fvals']);
end;

%% Report MCS results

K = results(:,1); mK = mean(K); delK = std(K);
lambdae = results(:,2); mle = mean(lambdae); delle = std(lambdae);
lambdab = results(:,3); mlB = mean(lambdab); dellB = std(lambdab);
ml = mle + mlB; dell = sqrt( delle.^2 + dellB.^2 );

% check this 
disp(['MCS K = ' sprintf('%0.4f',mK.*1e3) ' +/- ' sprintf('%0.4f',delK.*1e3) ' (' sprintf('%0.2f',100*(delK./mK)) '%)']);
disp(['MCS lambda_e = ' sprintf('%0.4f',mle.*1e10) ' +/- ' sprintf('%0.4f',delle.*1e10) ' (' sprintf('%0.2f',100*(delle./mle)) '%)']);
disp(['MCS lambda_B = ' sprintf('%0.4f',mlB.*1e10) ' +/- ' sprintf('%0.4f',dellB.*1e10) ' (' sprintf('%0.2f',100*(dellB./mlB)) '%)']);
disp(['MCS lambda = ' sprintf('%0.4f',ml.*1e10) ' +/- ' sprintf('%0.4f',dell.*1e10) ' (' sprintf('%0.2f',100*(dell./ml)) '%)']);

% report better age estimates for samples

for a = 2:(size(results,2)-3);
    thismt = mean(results(:,a+3));
    thisstd = std(results(:,a+3));
    disp([sprintf('%0.5f',thismt./1e6) ' +/- ' sprintf('%0.5f',thisstd./1e6) ' (' sprintf('%0.2f',100.*thisstd./thismt) ')']);
end;

% Figure and report covariances
cc = cov([K lambdae lambdab]);
disp('Covariances');
disp(['Cov (K,lambda_e) = ' num2str(cc(1,2))]);
disp(['Cov (K,lambda_b) = ' num2str(cc(1,3))]);
disp(['Cov (lambda_e,lambda_b) = ' num2str(cc(2,3))]);

%% 

%% plotting


% Set limits
minK = 1.6e-3; maxK = 1.7e-3;
minle = 5.6e-11; maxle = 6e-11;
minlB = 4.8e-10; maxlB = 5.1e-10;

% line plotting vectors
ks = linspace(minK,maxK,10);
les = linspace(minle,maxle,10);
lbs = linspace(minlB,maxlB,10);

% slice locations
bestK = xopt(1); bestle = xopt(2); bestlB = xopt(3); 


% --- K-lambda_e plot -------
figure(1);clf;
for a = 1:length(data.ti);
    % Set lB and plot K vs le
    yy = les;
    xx = (les./(data.ri(a).*(les+bestlB))).*(exp(data.ti(a).*(les+bestlB))-1);
    
    % error propagation
    dxxdti = (les./data.ri(a)).*(exp(data.ti(a).*(les+bestlB)));
    dxxdRi = (-les./((data.ri(a).^2).*(les+bestlB))).*(exp(data.ti(a).*(les+bestlB))-1);
    delxx = sqrt( (data.delri(a).*dxxdRi).^2 + (data.delti(a).*dxxdti).^2 );

    % plot error bounds
    xxx = [xx+delxx fliplr(xx-delxx)];
    yyy = [yy fliplr(yy)];
    patch(xxx,yyy,'g','tag','greenerror');hold on;
    % plot center line
    plot(xx,yy,'c','tag','greenline');
end;

% plot K, le determinations
plot([data.a(1) data.a(1)],[minle maxle],'b','tag','redline');
patch([data.a(1)-data.dela(1) data.a(1)-data.dela(1) data.a(1)+data.dela(1) data.a(1)+data.dela(1) data.a(1)-data.dela(1)],...
    [minle maxle maxle minle minle],'r','tag','rederror');
plot([minK maxK],[data.a(2) data.a(2)],'b','tag','redline');
patch([minK maxK maxK minK minK],...
    [data.a(2)-data.dela(2) data.a(2)-data.dela(2) data.a(2)+data.dela(2) data.a(2)+data.dela(2) data.a(2)-data.dela(2)],'r','tag','rederror');

% Data
plot(K,lambdae,'b.','tag','bluedot'); 
xlabel('K'); ylabel('lambda e');
plot(xopt(1),xopt(2),'ko','markerfacecolor','k','tag','blackdot');
axis([minK maxK minle maxle])

% ------ K - lambda B plot -----

figure(2);clf;
for a = 1:length(data.ti);
    % Set le and plot K vs lB
    yy = lbs;
    xx = (bestle./(data.ri(a).*(bestle+lbs))).*(exp(data.ti(a).*(bestle+lbs))-1);
    
    % error propagation
    dxxdti = (bestle./data.ri(a)).*(exp(data.ti(a).*(bestle+lbs)));
    dxxdRi = (-bestle./((data.ri(a).^2).*(bestle+lbs))).*(exp(data.ti(a).*(bestle+lbs))-1);
    delxx = sqrt( (data.delri(a).*dxxdRi).^2 + (data.delti(a).*dxxdti).^2 );
    
    % plot error bounds
    xxx = [xx+delxx fliplr(xx-delxx)];
    yyy = [yy fliplr(yy)];
    patch(xxx,yyy,'g','tag','greenerror');hold on;
    % plot center line
    plot(xx,yy,'c','tag','greenline');
end;

% plot K, lB determinations
plot([data.a(1) data.a(1)],[minlB maxlB],'b','tag','redline');
patch([data.a(1)-data.dela(1) data.a(1)-data.dela(1) data.a(1)+data.dela(1) data.a(1)+data.dela(1) data.a(1)-data.dela(1)],...
    [minlB maxlB maxlB minlB minlB],'r','tag','rederror');
plot([minK maxK],[data.a(3) data.a(3)],'b','tag','redline');
patch([minK maxK maxK minK minK],...
    [data.a(3)-data.dela(3) data.a(3)-data.dela(3) data.a(3)+data.dela(3) data.a(3)+data.dela(3) data.a(3)-data.dela(3)],'r','tag','rederror');

% Data
plot(K,lambdab,'b.','tag','bluedot'); 
xlabel('K'); ylabel('lambda B'); 
plot(xopt(1),xopt(3),'ko','markerfacecolor','k','tag','blackdot');
axis([minK maxK minlB maxlB])


% ------ le-lb plot -----------

figure(3);  clf; 
for a = 1:length(data.ti);
    % Set K and plot le vs lb
    xx = les;
    yy = zeros(size(xx));
    for b = 1:length(xx);
        % Find value of lb at central point
        yy(b) = fzero(@(x) (les(b)./(bestK.*(les(b)+x))).*(exp(data.ti(a).*(les(b)+x))-1) - data.ri(a),bestlB);
        % Find value of lb at 1 SD increase in ri
        upR(b) = fzero(@(x) (les(b)./(bestK.*(les(b)+x))).*(exp(data.ti(a).*(les(b)+x))-1) - (data.ri(a)+data.delri(a)),bestlB);
        % Find value of lb at 1 SD increase in ti
        upt(b) = fzero(@(x) (les(b)./(bestK.*(les(b)+x))).*(exp((data.ti(a)+data.delti(a)).*(les(b)+x))-1) - data.ri(a),bestlB);
    end;
    
    % Error propagation
    dyydRi = (upR-yy)./data.delri(a);
    dyydti = (upt-yy)./data.delti(a);
    delyy = sqrt( (dyydRi.*data.delri(a)).^2 + (dyydti.*data.delti(a)).^2 );
    
    % plot error bounds
    xxx = [xx fliplr(xx)];
    yyy = [yy+delyy fliplr(yy-delyy)];
    patch(xxx,yyy,'g','tag','greenerror');hold on;
    % plot center line
    plot(xx,yy,'c','tag','greenline');
    
end;

% plot le, lB determinations
plot([data.a(2) data.a(2)],[minlB maxlB],'b','tag','redline');
patch([data.a(2)-data.dela(2) data.a(2)-data.dela(2) data.a(2)+data.dela(2) data.a(2)+data.dela(2) data.a(2)-data.dela(2)],...
    [minlB maxlB maxlB minlB minlB],'r','tag','rederror');
plot([minle maxle],[data.a(3) data.a(3)],'b','tag','redline');
patch([minle maxle maxle minle minle],...
    [data.a(3)-data.dela(3) data.a(3)-data.dela(3) data.a(3)+data.dela(3) data.a(3)+data.dela(3) data.a(3)-data.dela(3)],'r','tag','rederror');

% plot ltot determination

xx = les;
yy = data.a(4)-xx;
xxx = [les fliplr(les)];
yyy = [yy+data.dela(4) fliplr(yy-data.dela(4))];
patch(xxx,yyy,[0.5 0.5 0.5],'tag','rederror');
plot(xx,yy,'k','tag','redline');


% data
plot(lambdae,lambdab,'b.','tag','bluedot'); 

xlabel('lambda e'); ylabel('lambda B'); 
plot(xopt(2),xopt(3),'ko','markerfacecolor','k','tag','blackdot');
axis([minle maxle minlB maxlB]);
  

%% Make bias plot

tU = data.ti;
deltU = data.delti;
%le = xopt(2); lb = xopt(3); K = xopt(1); 
% use starting estimates of parameters instead
% le = data.a(2);lb = data.a(3); K = data.a(1); 


tAr = (1./ml).*log(mK.*data.ri.*(ml./mle) + 1);
dtArdR = (1./ml).*(1./(data.ri.*mK.*(ml./mle) + 1)).*(mK.*ml./mle);
deltAr = sqrt( (dtArdR.*data.delri).^2);

bias = tU./tAr;
dbdtU = 1./tAr;
dbdtAr = tU./(tAr.^2);
dbias = sqrt( (dbdtU.*deltU).^2 + (dbdtAr.*deltAr).^2 );

figure(4); 
semilogx(tAr,bias,'ro');hold on;
for a = 1:length(tU);
    xx = [tAr(a) tAr(a)];
    yy = [bias(a)-dbias(a) bias(a)+dbias(a)];
    plot(xx,yy,'r');
end;
axis([1e7 1e10 0.99 1.03]);

%% Make normalized residuals plot

edges = -3:0.25:3;
ns = histc(misses,edges);
figure;
stairs(edges,ns,'g');
axis([-3 3 0 8]);
    
xn = -3:0.05:3;
y = 0.25.*length(misses).*exp(-0.5 * xn .^2) ./ (sqrt(2*pi));
hold on;
plot(xn,y,'r');




%% write out Monte Carlo results for main parameters


if saveFlag;
    fid = fopen([save_fname '.dat'],'w');
    fprintf(fid,'%s\n','K lambda_e lambda_b');
    fprintf(fid,'%0.5e %0.5e %0.5e\n',results(:,1:3)');
    fclose(fid);
end;
    
%% Make marginalized probability distrbution plots for three main
% parameters

K = results(:,1); mK = mean(K); delK = std(K);
lambdae = results(:,2); mle = mean(lambdae); delle = std(lambdae);
lambdab = results(:,3); mlB = mean(lambdab); dellB = std(lambdab);

minK = 1.61e-3; maxK = 1.67e-3;
minle = 5.65e-11; maxle = 5.85e-11;
minlB = 4.93e-10; maxlB = 5.02e-10;

nums = 5;
minK = mK-nums*delK; maxK = mK+nums.*delK;
minle = mle-nums*delle; maxle = mle+nums.*delle;
minlB = mlB-nums*dellB; maxlB = mlB+nums.*dellB;


numbins = 40;

figure; subplot(3,1,1);

edges1 = linspace(minK,maxK,numbins);
ns1 = histc(K,edges1);
stairs(edges1,ns1,'g');
set(gca,'xlim',[minK maxK],'ylim',[0 600]);
hold on;
plot([mK mK],[0 600],'k');
plot([mK-delK mK-delK mK+delK mK+delK mK-delK],[0 600 600 0 0],'b');

subplot(3,1,2);

edges2 = linspace(minle,maxle,numbins);
ns2 = histc(lambdae,edges2);
stairs(edges2,ns2,'g');
set(gca,'xlim',[minle maxle],'ylim',[0 600]);
hold on;
plot([mle mle],[0 600],'k');
plot([mle-delle mle-delle mle+delle mle+delle mle-delle],[0 600 600 0 0],'b');


subplot(3,1,3);

edges3 = linspace(minlB,maxlB,numbins);
ns3 = histc(lambdab,edges3);
stairs(edges3,ns3,'g');
set(gca,'xlim',[minlB maxlB],'ylim',[0 600]);
hold on;
plot([mlB mlB],[0 600],'k');
plot([mlB-dellB mlB-dellB mlB+dellB mlB+dellB mlB-dellB],[0 600 600 0 0],'b');



    



    
