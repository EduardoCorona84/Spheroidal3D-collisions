% Test multiparticle codes

%% Figure 1: Stokes suspension plot
% Generate large system of ~90 spheroid suspensions inside a large oblate
% shell, and simulate curvature flow using Stokes LP from Laplace LP.
% Currently plotting at set Z slices.

[Stk_XX,Stk_YY,Stk_ZZ]=test_stokes(100,8,0);

%% Figure 2,3
plot_sphpro_dirC();
plot_neuS_sphmixed(); % Boxchart neuS for prolates; 3D density plot for mixed system.

%% Figure 4: Condition number of completion to I/2 + D + C
f = figure;
f.Position=[100 300 1800 400];
subplot(1,3,1);
[errC,condC]=compare_aspect(0,0,0,0);

subplot(1,3,2);
[etaC,etaS]=scaled_C();

subplot(1,3,3);
scaled_S();

%% Solving BVP using BIE for Dirichlet or Neumann problems
boxchart_spectral(0,0,1); % C, dir, mixed
% boxchart_spectral(1,0,1); % S, dir, mixed
% boxchart_spectral(0,1,1); % C, neu, mixed
% boxchart_spectral(1,1,1); % S, neu, mixed
% boxchart_spectral(1,1,0);

% %% How condition number changes as aspect ratio grows
% f=figure;
% f.Position=[100 300 1400 400];
% subplot(1,3,1);
% [errC,condC]=compare_aspect(0,0,0,1);
% subplot(1,3,2);
% [errS,condS]=compare_aspect(1,0,0,1);
% subplot(1,3,3);
% [errSeta,condSeta]=compare_aspect(1,0,1,1);

%% Oblate eta scaling
[etac,etas]=scaled_C(1);
R_change=[1.01,1.1:0.1:4.1,5,6,7,8,10,100,500,1000];
% plot(R_change,etas); hold on;
% plot(R_change,R_change./4+1/4);
figure; plot(R_change,etac); 



%% Functions

function plot_sphpro_dirC()
    % First plot system of three prolates on which the BCP is solved for
    % p=24, distance = 0.1;
    % next plot boxchart_spectral(0,0) to get dirichlet BCP error with C as
    % completion term, and plot them in the same figure.

    % Figure 2: system of prolates
    parr = [4,8,12,16,24];  %p-values to test
    distance=1.*10.^(-(1:6));

    rng(13.5);

    f=figure;
    f.Position=[100 300 1150 400];
    subplot('Position',[0.05 .15 .3 .8])
    hold on;
    for l=1:length(parr)
        p=parr(l);
        np=2*p*(p+1);
        for j=1:length(distance)
            ifplot = (j==1)&&(l==length(parr));
            test_three_plot(p,10,[1.1 1.2 1.3],ifplot); 
        end
        grid on;
    end
   
    % Figure 3: Boxchart
    distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    errs= cell(1,length(parr));

    % write to file
    boxplotName='./spheroidal/data/boxplotData_useS0_ifNeu0.txt';

    fID = fopen(boxplotName,'r');

    for l=1:length(parr)
        line1=fgetl(fID);
        line1split=split(line1,'=');
        p=str2double(line1split(end));

        num_entries=str2double(fgetl(fID));
        err_p_1d=zeros(num_entries,1);

        for ind=1:num_entries
            err_p_1d(ind)=str2double(fgetl(fID));
        end

        err_p=reshape(err_p_1d,[],length(distance));
        errs{l} = err_p;

        fgetl(fID);
    end

    fclose(fID);

    for l=1:length(parr)
        % subplot(1,length(parr)+1,l)
        subplot('Position',[0.4+.12*(l-1) .15 .1 .8])
        hold on;
        ylim([-14,2])
        boxchart(errs{l});
        xlabel("-log10(distance)")
        if l==1
            ylabel("log10 relative error")
        end
        title("p="+string(parr(l)));
        fontsize(gca,13,'points');
        hold off;
    end

end

function plot_neuS_sphmixed()

    f=figure;
    f.Position=[100 300 1150 400];

    subplot('Position',[0.05 .15 .3 .8])
    parr = [4,8,12,16,24];  %p-values to test
    distance=1.*10.^(-(1:6));
    rng(13.5);
    hold on;
    for l=1:length(parr)
        p=parr(l);
        np=2*p*(p+1);
        for j=1:length(distance)
            ifplot = (j==1)&&(l==length(parr));
            test_three_plot(p,10,[1.1 1.2 1.3],ifplot,1); 
        end
        grid on;
    end

    % Boxchart for exterior Neuman with S, mixed system
    parr = [4,8,12,16,24];  %p-values to test
    distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    errs= cell(1,length(parr));

    % recover from file
    boxplotName='./spheroidal/data/boxplotData_useS1_ifNeu1_mixed.txt';

    fID = fopen(boxplotName,'r');

    for l=1:length(parr)
        line1=fgetl(fID);
        line1split=split(line1,'=');
        p=str2double(line1split(end));

        num_entries=str2double(fgetl(fID));
        err_p_1d=zeros(num_entries,1);

        for ind=1:num_entries
            err_p_1d(ind)=str2double(fgetl(fID));
        end

        err_p=reshape(err_p_1d,[],length(distance));
        errs{l} = err_p;

        fgetl(fID);
    end

    fclose(fID);

    for l=1:length(parr)
        % subplot(1,length(parr)+1,l)
        subplot('Position',[0.4+.12*(l-1) .15 .1 .8])
        hold on;
        ylim([-14,2])
        boxchart(errs{l});
        xlabel("-log10(distance)")
        if l==1
            ylabel("log10 relative error")
        end
        title("p="+string(parr(l)));
        fontsize(gca,13,'points');
        hold off;
    end



    % 
    % % Boxchart for exterior Dirichlet with S as completion term
    % parr = [4,8,12,16,24];  %p-values to test
    % distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    % errs= cell(1,length(parr));

    % recover from file
    % boxplotName='./spheroidal/data/boxplotData_useS1_ifNeu1_mixed.txt';
    % 
    % fID = fopen(boxplotName,'r');
    % 
    % for l=1:length(parr)
    %     line1=fgetl(fID);
    %     line1split=split(line1,'=');
    %     p=str2double(line1split(end));
    % 
    %     num_entries=str2double(fgetl(fID));
    %     err_p_1d=zeros(num_entries,1);
    % 
    %     for ind=1:num_entries
    %         err_p_1d(ind)=str2double(fgetl(fID));
    %     end
    % 
    %     err_p=reshape(err_p_1d,[],length(distance));
    %     errs{l} = err_p;
    % 
    %     fgetl(fID);
    % end
    % 
    % fclose(fID);
    % 
    % for l=1:length(parr)
    %     % subplot(1,length(parr)+1,l)
    %     subplot('Position',[0.05+.12*(l-1) .15 .1 .8])
    %     hold on;
    %     ylim([-18,2])
    %     boxchart(errs{l});
    %     xlabel("-log10(distance)")
    %     if l==1
    %         ylabel("log10 relative error")
    %     end
    %     title("p="+string(parr(l)));
    %     fontsize(gca,13,'points');
    %     hold off;
    % end

    % % System of 2 oblates, 1 prolate
    % parr = [4,8,12,16,24];  %p-values to test
    % distance=1.*10.^(-(1:6));
    % 
    % rng(13.5);
    % subplot('Position',[.68 .15 .3 .8])
    % hold on;
    % for l=1:length(parr)
    %     p=parr(l);
    %     np=2*p*(p+1);
    %     for j=1:length(distance)
    %         ifplot = (j==1)&&(l==length(parr));
    %         test_three_plot(p,10,[1.1 1.2 1.3],ifplot,1); 
    %     end
    %     grid on;
    % end

end

% Previously fig23()
function plot_dirS_aR()
    % Run same code as boxchart_spectral(0,1) to get neumann BCP error with C as
    % completion term, and run compare_aspect(0,0) to get effect of aspect
    % ratios, and plot them in the same figure.

    % Figure 5: Boxchart
    parr = [4,8,12,16,24];  %p-values to test
    distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    errs= cell(1,length(parr));

    % write to file
    boxplotName='./spheroidal/data/boxplotData_useS0_ifNeu1.txt';

    fID = fopen(boxplotName,'r');

    for l=1:length(parr)
        line1=fgetl(fID);
        line1split=split(line1,'=');
        p=str2double(line1split(end));

        num_entries=str2double(fgetl(fID));
        err_p_1d=zeros(num_entries,1);

        for ind=1:num_entries
            err_p_1d(ind)=str2double(fgetl(fID));
        end

        err_p=reshape(err_p_1d,[],length(distance));
        errs{l} = err_p;

        fgetl(fID);
    end

    fclose(fID);
    
    f=figure;
    f.Position=[100 300 1150 400];
    for l=1:length(parr)
        % subplot(1,length(parr)+1,l)
        subplot('Position',[0.05+.12*(l-1) .15 .1 .8])
        hold on;
        ylim([-14,2])
        boxchart(errs{l});
        xlabel("-log10(distance)")
        if l==1
            ylabel("log10 relative error")
        end
        title("p="+string(parr(l)));
        fontsize(gca,13,'points');
        hold off;
    end

    % Figure 6: Aspect ratio
    R_change=[1.1,1.5,2,4];
    errs_1 = cell(1,length(R_change));
    set_dist = 0.5;
    Xeval_p=8; np_eval=2*Xeval_p*(Xeval_p+1);

    rng(7.689);

    for l=1:length(parr)
        p=parr(l);
        this_p_errs=zeros(np_eval,length(R_change));
        for j=1:length(R_change)
            this_u0=1/sqrt(1-1/R_change(j)^2);
            u0_list=[this_u0];
            [soln,truesoln,~] = test(p,10,1,u0_list,set_dist,0,0,0,1,Xeval_p);
            this_p_errs(:,j) = log10(abs((soln-truesoln)./max(truesoln)));

        end
        errs_1{l}=this_p_errs;
    end
    
    % subplot(1,length(parr)+1,length(parr)+1)
    subplot('Position',[.68 .15 .3 .8])
    ylim([-22,2]); 
    hold on;
    legends=[];
    for l=1:length(parr)
        boxchart(errs_1{l});
        avg_col=mean(errs_1{l},1);
        plot(avg_col,'-o',"MarkerSize",10,'lineWidth',1.25);
        legends=[legends 'p='+string(parr(l)) 'Mean p='+string(parr(l))];
    end
    xval=cell(1,length(R_change));
    log2R=round(log2(R_change),2);
    for j=1:length(R_change)
        xval{j}=compose("%.2f",log2R(j));
    end
    yval=cell(1,13);
    for j=1:13
        yval{j}=string(-22+2*(j-1));
    end
    set(gca,'XTickLabel',xval);
    set(gca,'YTickLabel',yval);
    set(gca,'YTick', -22:2:2);
    xlabel("log2(aspect ratio)")
    ylabel("log10 relative error")
    fontsize(gca,13,'points');
    lgd = legend(legends);
    fontsize(lgd,12,'points');
    set(lgd,'position',[.91 .27 .05 .15]); % also need to drag away horizontal white space.
    hold off;
end

function errs = boxchart_spectral(useS,if_neumann,mix_obl)
    % Set up: 3 prolate spheroids in space, with 2 random point charges
    % inside each. Use Dirichlet or Neumann problem formulation.
    % Demonstrates spectral convergence and maintained accuracy for
    % near-singular evalulation.
    if nargin==2
        mix_obl = 0;
    end

    ns=3;
    % parr = [4,8,12,16];  %p-values to test
    parr = [4,8,12,16,24];
    distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    errs= cell(1,length(parr));

    % write to file
    boxplotName='./spheroidal/data/boxplotData_useS'+string(useS)+'_ifNeu'+string(if_neumann)+'.txt';

    rFlag=exist(boxplotName,'file');
    % rFlag=0;
    if rFlag
        fID = fopen(boxplotName,'r');

        for l=1:length(parr)
            line1=fgetl(fID);
            line1split=split(line1,'=');
            p=str2double(line1split(end));

            num_entries=str2double(fgetl(fID));
            err_p_1d=zeros(num_entries,1);

            for ind=1:num_entries
                err_p_1d(ind)=str2double(fgetl(fID));
            end

            err_p=reshape(err_p_1d,[],length(distance));
            errs{l} = err_p;

            fgetl(fID);
        end

        fclose(fID);

    else
        rng(13.5);

        for l=1:length(parr)
            p=parr(l);
            np=2*p*(p+1);
            this_p_errs_spec=zeros(np*ns,length(distance));
            for j=1:length(distance)
                [soln,truesoln,~,~] = test(p,5,3,[1.1 1.2 1.3],distance(j).*ones(1,3),useS,0,if_neumann,1,p,mix_obl);
                this_p_errs_spec(:,j) = log10(abs((soln-truesoln)./max(truesoln)));
            end
            errs{l}=this_p_errs_spec;
        end

        fID=fopen(boxplotName,'w');
        for l=1:length(parr)
            fprintf(fID,'p=%d\n',parr(l));
            err_p_1d=reshape(errs{l},[],1);
            fprintf(fID,'%d\n',length(err_p_1d));
            fprintf(fID,'%f\r\n',err_p_1d);
            fprintf(fID,'\n');
        end
        fclose(fID);

    end
    
    % figure;
    % for l=1:length(parr)
    %     subplot(1,length(parr),l)
    %     hold on;
    %     ylim([-14,2])
    %     boxchart(errs{l});
    %     xlabel("-log10(distance)")
    %     if l==1
    %         ylabel("log10 relative error")
    %     end
    %     title("p="+string(parr(l)));
    %     fontsize(gca,15,'points');
    %     hold off;
    % end
end

% plot boxchart error and condition number (if indicated).
function [errs_1,cond_nums] = compare_aspect(useS,if_neumann,scaled, plt_cond,obl)
    % Compare accuracy of method when aspect ratio of the prolate spheroid
    % increases.
    % Set up: 1 prolate spheroid centered at the origin, with 2 random point charges
    % inside slightly off the z-axis. Use Dirichlet or Neumann problem formulation.

    if nargin==4
        obl=0;
    end

    parr = [4,8,12,16,24];  %p-values to test
    R_change=[1.1,1.5,2,4,8];
    % R_change=[1.1,1.5,2,4,8,20,50];
    
    errs_1 = cell(1,length(R_change));
    cond_nums=zeros(length(parr),length(R_change));

    set_dist = 0.5;
    Xeval_p=16; np_eval=2*Xeval_p*(Xeval_p+1);

    rng(7.689);

    % set up x axis 
    xval=cell(1,length(R_change));
    log2R=round(log2(R_change),2);
    for j=1:length(R_change)
        xval{j}=compose("%.2f",log2R(j));
    end

    % figure;
    for l=1:length(parr)
        p=parr(l);
        this_p_errs=zeros(np_eval,length(R_change));

        for j=1:length(R_change)

            if obl
                aR = R_change(j);
                this_u0 = aR./sqrt(1-aR.^2);
            else
                this_u0=1/sqrt(1-1/R_change(j)^2);
            end

            if scaled && useS
                eps_sbt=1/R_change(j);
                scale=1/(2*eps_sbt*log(1/eps_sbt));
            else
                scale=1;
            end

            [soln,truesoln,~,condK] = test(p,5,1,[this_u0],set_dist,useS,0,if_neumann,scale,Xeval_p,1);

            this_p_errs(:,j) = log10(abs((soln-truesoln)./max(truesoln)));

            cond_nums(l,j)=condK;

        end
        errs_1{l}=this_p_errs;

        if plt_cond
            plot(log2R,cond_nums(l,:),"-*"); hold on;
        end
    end

    if plt_cond
        xticks(log2R)
        set(gca,'XTickLabel',xval);
        set(gca,'TickLabelInterpreter','latex');
        xlabel("log2(aspect ratio)")
        ldg=legend("p=4","p=8","p=12","p=16");
        
        if useS
            if scaled
                set(gca,'TickLabelInterpreter','latex');
                ylabel("cond(1/2 I + D + \eta S)")
                set(ldg,'Position',[0.8,0.785,0.05,0.1]);
            else
                ylabel("cond(1/2 I + D + S)")
                set(ldg,'Position',[0.5,0.785,0.05,0.1]);
            end
        else
            ylabel("cond(1/2 I + D + C)")
            set(ldg,'Position',[0.2,0.785,0.05,0.1]);
        end
    end

    
    % Boxchart of errors at different p, as aspect ratio increases.
    % figure;
    % hold on;
    ylim([-22,2]);
    legends=[];
    for l=1:length(parr)
        boxchart(errs_1{l}); hold on;
        avg_col=mean(errs_1{l},1);
        plot(avg_col,'-o',"MarkerSize",10,'lineWidth',1.25);
        legends=[legends 'p='+string(parr(l)) 'Mean p='+string(parr(l))];
    end

    yval=cell(1,13);
    for j=1:13
        yval{j}=string(-22+2*(j-1));
    end
    % set(gca,'TickLabelInterpreter','latex');
    set(gca,'XTickLabel',xval);
    set(gca,'YTickLabel',yval);
    set(gca,'YTick', -22:2:2);
    xlabel("log2 aspect ratio")
    ylabel("log10 relative error")
    fontsize(gca,15,'points');
    lgd = legend(legends);
    fontsize(lgd,15,'points');
    set(lgd,'position',[0.3,0.27,0.05,0.1]) % also need to drag away horizontal white space.
    hold off;
end

% reports the argmin and min of cond(K) as function of scalar eta. 
function [etaClst,etaSlst]=scaled_C(obl)

    if nargin==0
        obl=0;
    end

    % parr = [4,8,12,16];  % cond() for all p values look the same
    parr = [8];
    if ~obl
        R_change=[1.1,1.5,2,4,8,16,32,64];
    else
        R_change=[1.01,1.1:0.1:4.1,5,6,7,8,10,100,500,1000];
    end

    cond_C=zeros(length(parr),length(R_change));
    cond_etaC = cond_C;
    cond_S = cond_C;
    cond_etaS = cond_C;

    cond_etaC2=cond_C;

    etaClst = zeros(size(R_change));
    etaSlst = etaClst; etaC2lst=etaClst;

    % set up x axis 
    xval=cell(1,length(R_change));
    log2R=round(log2(R_change),2);
    for j=1:length(R_change)
        xval{j}=compose("%.2f",log2R(j));
    end

    figure;
    for l=1:length(parr)
        p=parr(l); np = 2*p*(p+1);

        for j=1:length(R_change)

            if ~obl
                [CM,SM,DM,~] = get_mats(1/R_change(j),p);
            else
                [CM,SM,DM,~] = get_mats_ob(1/R_change(j),p);
            end

            options = optimoptions('fmincon','Display','off');
            [etaC,condKC] = fmincon(@(s) cond(.5*eye(np)+DM+s*CM),1,-1,-1e-10,[],[],[],[],[],options);
            [etaS,condKS] = fmincon(@(s) cond(.5*eye(np)+DM+s*SM),1,-1,-1e-10,[],[],[],[],[],options);

            cond_C(l,j) = cond(.5*eye(np)+DM+CM);
            cond_etaC(l,j)=condKC;
            cond_S(l,j) = cond(.5*eye(np)+DM+SM);
            cond_etaS(l,j) = condKS;

            ep=1/R_change(j);
            E=sqrt(1-ep.^2); 
            if ~obl
                SA = 2*pi*(ep.^2).*(1+(R_change(j)./E).*asin(E));
            else
                SA = 2*pi*(ep.^2).*(1+R_change(j)^2./2./E.*log((1+E)./(1-E)));
            end
            
            if R_change(j)<=3
                etaC2 = ep./SA; 
                % fprintf("\n etaC2 with SA=%.3g", etaC2);
            else
                etaC2=0.06;
            end
            etaC2lst(j)=etaC2;
            cond_etaC2(l,j) = cond(.5*eye(np)+DM+etaC2*CM);


            etaClst(j) = etaC;
            etaSlst(j) = etaS;

        end

        semilogy(log2R,cond_C(l,:),"-*","MarkerSize",10,'lineWidth',1.25); hold on;
        semilogy(log2R,cond_etaC(l,:),"--","MarkerSize",10,'lineWidth',1.25); hold on;
        semilogy(log2R,cond_S(l,:),"-o","MarkerSize",10,'lineWidth',1.25);
        semilogy(log2R,cond_etaS(l,:),"-.","MarkerSize",10,'lineWidth',1.25);
        % semilogy(log2R,cond_etaC2(l,:),"-s","MarkerSize",10,'lineWidth',1.25); 

        % semilogy(1:numel(log2R),cond_C(l,:),"-*","MarkerSize",10,'lineWidth',1.25); hold on;
        % semilogy(1:numel(log2R),cond_etaC(l,:),"--","MarkerSize",10,'lineWidth',1.25); hold on;
        % semilogy(1:numel(log2R),cond_S(l,:),"-o","MarkerSize",10,'lineWidth',1.25);
        % semilogy(1:numel(log2R),cond_etaS(l,:),"-.","MarkerSize",10,'lineWidth',1.25);

    end
    % xlim([0.5,7.5]);
    set(gca,'XTickLabel',xval);
    % set(gca,'TickLabelInterpreter','latex');
    xlabel("log2 aspect ratio ")
    lgd=legend("$\mathcal{C} = \mathcal{C}_I$, p=16","$\mathcal{C} = \eta \mathcal{C}_I$","$\mathcal{C} = \mathcal{S}$","$\mathcal{C} = \eta \mathcal{S}$","Interpreter",'latex');
    % lgd=legend("$\mathcal{C} = \mathcal{C}_I$, p=16","$\mathcal{C} = \eta \mathcal{C}_I$","$\mathcal{C} = \mathcal{S}$","$\mathcal{C} = \eta \mathcal{S}$","$\mathcal{C} = \varepsilon/S \ \mathcal{C}_I$","Interpreter",'latex');
    % lgd=legend("$\mathcal{CP} = \eta \mathcal{C}_I$","$\mathcal{CP} = \varepsilon/S \ \mathcal{C}_I$","Interpreter",'latex');
    ylabel("cond(I/2 + D + CP)")
    fontsize(gca,15,'points');
    fontsize(lgd,15,'points');
    set(lgd,'position',[0.435,0.775,0.05,0.1])
    hold off;

    keyboard;

end

% plots argmin for cond(K) and compare to CSBT value.
function scaled_S()
    parr = [16];
    R_change = [1.1:0.02:1.4,1.5,2,4,8,16,32,64];
    x_list=[1.1,1.5,2,4,8,16,32,64];

    etaSlst=zeros(size(R_change));

    % set up x axis 
    xval=cell(1,length(x_list));
    log2R=round(log2(x_list),2);
    for j=1:length(x_list)
        xval{j}=compose("%.2f",log2R(j));
    end

    % figure;
    for l=1:length(parr)
        p=parr(l); np = 2*p*(p+1);

        for j=1:length(R_change)

            [CM,SM,DM,~] = get_mats(1/R_change(j),p);
            % [CM,SM,DM,~] = get_mats_ob(1/R_change(j),p);

            options = optimoptions('fmincon','Display','off');
            [etaS,~] = fmincon(@(s) cond(.5*eye(np)+DM+s*SM),1,-1,-1e-10,[],[],[],[],[],options);

            etaSlst(j) = etaS;

        end

    end
    eps_sbt=1./R_change;
    scale=1/2./(eps_sbt.*log(1./eps_sbt)); 

    plot(log2(R_change),etaSlst,"-*","MarkerSize",10,'lineWidth',1.25); hold on;
    plot(log2(R_change),scale,'-o',"MarkerSize",10,'lineWidth',1.25); hold off;
    % xlim([0,6.5]);
    set(gca,'XTick',log2R);
    set(gca,'XTickLabel',xval);
    xlabel("log2 aspect ratio");
    lgd = legend("minimizer S","CSBT");
    fontsize(gca,15,'points');
    fontsize(lgd,15,'points');
    set(lgd,'Position',[0.715,0.775,0.05,0.1]);
    ylabel('$\eta$','Interpreter','latex');

end
    

function [soln,truesoln,sigma_vec,condK] = test(p,eta,ns,u0,target_distances,useS,plt,neumann,S_scale,Xeval_p,mix_obl)
    % (1) set up system of 3 spheroids, 2 close and 1 far
    % (2) put a few point charges around the center of each spheroid.
    %       - Calculate the potential explicitly.
    %       - Determine potential on the boundary of each spheroid for BCs, f
    % (3) Use MatVec self-evaluation (no X) to construct matrix operator D.
    % (4) Solve (1/2*I + D + Completion)*sigma = f
    %       * Completion term will need some care, since it's multiple particles.
    %       For each particle it will be something like sum_i[ 1/|| x- c_i||
    %       integral(sigma_i) ]
    % (5) Compare sigma with true point charge potential at targets that
    %     are 'target_distances' away from surface of each spheroid.

    if nargin==8
        S_scale=1;
    end
    if nargin==9
        Xeval_p=p;
    end
    if nargin==10
        mix_obl = 0;
    end

    % Set up spheroid system
    np=2*p*(p+1);
    
    if length(u0)~=ns
        fprintf("\n input u0 length does not match number of spheroids indicated.");
        ns=length(u0);
    end
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal=0;
    pars.u0=u0;

    if ns==3
        if mix_obl
            obl = [1,0,1];
        else
            obl = [0,0,0];
        end
        pars.oblate = obl;
        a = 1./u0;
        a(obl==1) = 1./sqrt(u0(obl==1).^2+1);
        pars.a = a;
        
        pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
        % pars.centers = [0 0 0; 2.5 0 0; 1.6 1.6 1.6];
        thetas = [0 pi/10 5*pi/3];
        phis = [0 0 pi/5];
    else
        if ns~=1
            fprintf("\n Number of spheroids given not implemented here. Setting ns=1.")
            ns=1;
        end
        if mix_obl
            pars.a=1/sqrt(u0^2+1);
            pars.oblate=1;
        else
            pars.a=1/u0;
            pars.oblate=0;
        end
        
        pars.centers=[0 0 0];
        thetas=0;
        phis=0;
    end

    Ri=zeros(3,3,ns);
    for i=1:ns
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,ns);

    % pars.plot();
        
    % Point charges
    % rng(13.5) 
    nc = 2;
    c=pars.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % random point charges
    if ns==3
        scale = repmat(.5.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
        scale = scale(:);
        d = repmat(scale,1,3).*(rand(size(Xptch))-.5);
    else
        d_z=0.5.*(rand(nc,1)-.5);
        d_xy=0.01.*(rand(nc,2)-0.5);
        d = [d_xy d_z];
    end
    Xptch = Xptch + d; %point charge locations, near centers of spheroids
    ptch = (2.*rand(ns*nc,1)-1); %charge value, between -1 and 1

    % % determined point charges
    % ptch = [0.1;0.6];
    
    if plt
        pars.plot(Xptch);
        title("spheroids and point charges")
    end

    Y=pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p);
    
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    
    if ~neumann
        if useS
            CM = S_scale.*spheroidalMatVecKernel(pars,'SL',p);
        else
            CM=[];
            for i=1:ns
                if ~pars.oblate(i)
                    Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                else
                    Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                end
                ci=c(i,:);
                
                % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
                % respective spheroid surface
                Sns = SurfaceSph(Yi(:));
                W = Sns.geoProp.W;
            
                RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
                CM = [CM (1./RY)*(W'.*wt)];
            end
        end
        % We will use this for surface boundary conditions
        truesolnSurf=PtChargePotential(ptch,Xptch,Y);

        % Construct DL on-surface matrices
        DM = spheroidalMatVecKernel(pars,'DL',p);
        
        % Final operator: 1/2 I + DM + CM
        K = .5*eye(ns*np) + DM +CM;
    else
        truesolnSurf=PtChargeFlux(ptch, Xptch,Y,NrY);

        % Construct SP on-surface matrices
        SPM= spheroidalMatVecKernel(pars,'SP',p);
        
        % Final operator: -1/2 I + SPM
        K = -.5*eye(ns*np) + SPM;
    end

    condK=cond(K);
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
    
    % Target points to test at
    [theta,phi]=gl_grid(Xeval_p);
    v=cos(theta);
    Xcell=cell(1,ns);
    for i=1:ns
        u0=pars.u0(i);
        if ~pars.oblate(i)
            normal=1./sqrt((u0^2-1).*(u0^2-v.^2)).*[u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*sin(phi) (u0.^2-1).*v];
            Yi = prolate_spheroid_shape(Xeval_p,pars.u0(i),pars.a(i),'cart');
        else
            normal=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
            Yi = oblate_spheroid_shape(Xeval_p,pars.u0(i),pars.a(i),'cart');
        end
        Xcell{i} = Yi + target_distances(i).*normal;
    end
    Xeval = pars.set_X_targets(Xcell);

    if plt
        pars.plot(Xeval);
        title('spheroids and target points')
    end
    
    if ~neumann
        if useS
            C = S_scale.*spheroidalMatVec(pars,'SL',Xeval);
        else
            C=zeros(size(Xeval,1),1);
            for i=1:ns
            
                if ~pars.oblate(i)
                    Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                else
                    Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                end
                ci=c(i,:);
                
                % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
                % respective spheroid surface
                Sns = SurfaceSph(Yi(:));
                sigmaSurfInt = integrateOverS(Sns,sigma(:,:,i));
            
                RXeval = vecnorm(Xeval-repmat(ci,size(Xeval,1),1),2,2);
                
                C = C+ sigmaSurfInt./RXeval;
            end
        end

        soln=spheroidalMatVec(pars,'DL',Xeval) + C;

    else

        soln=spheroidalMatVec(pars,'SL',Xeval);

    end

    truesoln=PtChargePotential(ptch,Xptch,Xeval); 

    fprintf('p=%d: spectral method err=%.6e\n',p, max(abs(truesoln-soln)))  

end



function test_three_plot(p,eta,u0,plt,mix_obl)
    if nargin==4
        mix_obl = 0;
    end

    % Set up spheroid system
    np=2*p*(p+1);
    ns=3;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal=0;
    pars.u0=u0;

    if mix_obl
        obl = [1,0,1];
    else
        obl = [0,0,0];
    end
    pars.oblate = obl;
    a = 1./u0;
    a(obl==1) = 1./sqrt(u0(obl==1).^2+1);
    pars.a = a;

    % pars.a = 1./u0;
    % pars.oblate = [0 0 0];
    
    pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
    thetas = [0 pi/10 5*pi/3];
    phis = [0 0 pi/5];
    Ri=zeros(3,3,3);
    for i=1:3
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,ns);
        
    % Point charges
    nc = 2;
    scale = repmat(.5.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
    scale = scale(:);
    c=pars.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % random point charges
    d = repmat(scale,1,3).*(rand(size(Xptch))-.5);
    Xptch = Xptch + d; %point charge locations, near centers of spheroids
    ptch = (2.*rand(ns*nc,1)-1); %charge value, between -1 and 1

    if ~plt
        return;
    end

    Y=pars.get_X;

    % truesolnSurf=PtChargePotential(ptch,Xptch,Y);

    Y=pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p);
    
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    CM=[];
    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        ci=c(i,:);
        
        % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
        % respective spheroid surface
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;
    
        RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
        CM = [CM (1./RY)*(W'.*wt)];
    end

    % We will use this for surface boundary conditions
    truesolnSurf=PtChargePotential(ptch,Xptch,Y);

    % Construct DL on-surface matrices
    DM = spheroidalMatVecKernel(pars,'DL',p);
    
    % Final operator: 1/2 I + DM + CM
    K = .5*eye(ns*np) + DM +CM;
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;



    addpath('spheroidal/plot/utils/');
    cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/python3"');

    for sphind=1:3
        Yi=Y((sphind-1)*np+1:sphind*np,:);
        % shades=truesolnSurf((sphind-1)*np+1:sphind*np,:);
        shades=pars.sigma(:,:,sphind);
        plotb([Yi(:,1);Yi(:,2);Yi(:,3)],real(shades)); hold on;
        colormap(cmp)
    end
    view(45,10);
    cb=colorbar; 
    % if mix_obl
    %     set(cb,'Position',[0.95 0.3 0.01 0.3]);
    % else
        set(cb,'Position',[0.32 0.30 0.01 0.3]);
    % end
    set(cb,'YTick',-0.1:0.1:0.1);
    hold off;
end


function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
% np=length(Y);
np=size(Y,1);
Rptch=zeros(np,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end

function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

function [CM,SM,DM,this_u0] = get_mats(eps,p)

    np=2*p*(p+1);

    R = 1/(eps);

    pars=SpheroidalParameters;

    this_u0=1/sqrt(1-1/R^2);
    pars.u0=this_u0;
    % pars.a=1/this_u0;
    pars.a = 20/this_u0;
    pars.oblate=0;

    pars.matvec_eta = 5;
    pars.isReal=0;
    
    pars.centers=[0 0 0];
    thetas=0;
    phis=0;
    Ri=zeros(3,3,1);
    for i=1:1
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,1); pars.get_shc();
    
    SM = spheroidalMatVecKernel(pars,'SL',p);

    DM = spheroidalMatVecKernel(pars,'DL',p);

    CM=[];
    Yi = prolate_spheroid_shape(p,pars.u0,pars.a,'cart');
    
    % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
    % respective spheroid surface
    Sns = SurfaceSph(Yi(:));
    W = Sns.geoProp.W;

    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    RY = vecnorm(Yi,2,2);
    CM = [CM (1./RY)*(W'.*wt)];

end

function [CM,SM,DM,this_u0] = get_mats_ob(eps,p)
    np=2*p*(p+1);

    R = 1/(eps);

    pars=SpheroidalParameters;

    % this_u0=1/sqrt(1-1/R^2);
    % pars.u0=this_u0;
    % pars.a=1/this_u0;

    this_u0 = 1/sqrt(R^2-1);
    pars.u0 = this_u0;
    pars.a = 1/sqrt(this_u0^2+1);

    pars.matvec_eta = 5;
    pars.isReal=0;
    pars.oblate=1;
    pars.centers=[0 0 0];
    thetas=0;
    phis=0;
    Ri=zeros(3,3,1);
    for i=1:1
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,1); pars.get_shc();
    
    SM = spheroidalMatVecKernel(pars,'SL',p);

    DM = spheroidalMatVecKernel(pars,'DL',p);

    CM=[];
    Yi = oblate_spheroid_shape(p,pars.u0,pars.a,'cart');
    
    % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
    % respective spheroid surface
    Sns = SurfaceSph(Yi(:));
    W = Sns.geoProp.W;

    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    RY = vecnorm(Yi,2,2);
    CM = [CM (1./RY)*(W'.*wt)];

end

