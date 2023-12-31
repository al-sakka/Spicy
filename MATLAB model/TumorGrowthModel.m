 clc
close all 
clearvars

%I define variables
G=4;                    % number of domains
n=zeros(201,201,G);   % total tumor cells
Ce=zeros(201,201,G);  % epithelial-like cells
Cm=zeros(201,201,G);  % mesenchymal-like cells
chosen=zeros(201,201,G);
agecm=zeros(201,201,G);   %Cells' age
agece=zeros(201,201,G);
dead=zeros(201,201,G);
deadage=zeros(201,201,G);
w=zeros(201,201,G);   % ECM density
neww=zeros(201,201,G);
m=zeros(201,201,G);   % MMP-2 concentration
newm1=zeros(201,201,G);
newm2=zeros(201,201,G);
blood=zeros(201,201,G);  % total Blood vessels
u=zeros(201,201,G);   % normal Blood vessels
v=zeros(201,201,G);   % reuptured Blood vessels
Vcm=zeros(201,201); % mesenchymal-like cells in vessels
Vce=zeros(201,201); % epithelial-like cells in vessels
agev=zeros(201,201);% age in the vessel
%II. Initial conditions
deltat=0.001;       % time step
deltax=0.005;       % space step
step=48000;
Dcm=0.0001;         % Mesenchymal-like cancer cell diffusion coefficient
Dce=0.00005;        % Epithelial-like cancer cell diffusion coefficient
phicm=0.0005;       % Mesenchymal haptotactic sensitivity coefficient
phice=0.0005;       % Epithelial haptotactic sensitivity coefficient
Dm=0.001;           % MMP-2 diffusion coefficient
theta=0.195;        % MMP-2 production rate
lambda=0.1;         % MMP-2 decay rate 
alpha1=1;           % ECM degradation rate by MT1-MMP
alpha2=1;           % ECM degradation rate by MMP-2
Tv=18;             % Time CTCs spend in the vasculature
Te=300;            % Epithelial doubling time
Tm=200;            % Mesenchymal doubling time
Td=0;            % time of recovery of the dead cells
Tdm=50000;           % life time of Mesenchymal cells
Tde=50000;           % life time of Epithelial cells
Ps=0.0005;          % Single CTC survival probability
Pc=0.025;           % CTC cluster survival probability
ebslon1=0.5461;     % Extravasation probability to first site
ebslon2=0.2553;     % Extravasation probability to second site
ebslon3=0.1986;     % Extravasation probability to third site
% ~~~~primary tumor~~~~~~%  
for x=1:201
    for y=1:201
            if x>=96 && x<=105 && y>=96 && y<=105
                n(x,y,1)=randi([0,4]);
                key=rand(1);
                if key<=0.4 && n(x,y,1)>0
                     Ce(x,y,1)=n(x,y,1);
                     agece(x,y,1)=floor(1+rand(1)*500);
                else 
                if key<=1 && n(x,y,1)>0
                     Cm(x,y,1)=n(x,y,1);
                     agecm(x,y,1)=floor(1+rand(1)*300);
                end
                end
            end
            w(x,y,:)=1;
    end
end
% ~~~~vessels~~~~~~%
for x=1:100
    for g=1:4
         if g==1
             i=randi(201);
             j=randi(201);
            if ((i>=2 && i<=94) || (i>=107 && i<=199))&& ((j>=2 && j<=94) || (j>=107 && j<=199)) && sum(u(:,:,g),"all")<8
                u(i,j,g)=randi([0,1]);
            end
            if ((i>=2 && i<=94) || (i>=107 && i<=199))&& ((j>=2 && j<=94) || (j>=107 && j<=199)) && sum(v(:,:,g),"all")<10 && u(i,j,g)==0
                p=randi([0,1]);
                v(i,j,g)=p;
                v((i+1),j,g)=p;
                v(i,(j+1),g)=p;
                v((i-1),j,g)=p;
                v(i,(j-1),g)=p;
            end
         else
            i=randi(201);
            j=randi(201);
            if ((i>=2 && i<=94) || (i>=107 && i<=199))&& ((j>=2 && j<=94) || (j>=107 && j<=199)) && sum(u(:,:,g),"all")<10
                u(i,j,g)=randi([0,1]);
            end
        end
    end
end
blood=u+v;
% ~~~~simulation~~~~~~%
for st=1:step
    st
    tempx=0;
    tempy=0;
    tempn=0;
    tempr=0;
    for x=1:201
    for y=1:201
    for g=1:G
        if n(x,y,g)>0
            tempx=tempx+x;
            tempy=tempy+y;
            tempn=tempn+n(x,y,g);
        end
        if dead(x,y,g)>0
            tempr=tempr+dead(x,y,g);
        end
    end
    end
    end
tempx=tempx/tempn;
tempy=tempy/tempn;
tempr=sqrt((tempn+tempr)/pi);
%1. MMP-2 All Sites
for g=1:G
    newm1(1,:,g)=(1-Dm)*m(1,:,g)+Dm*m(1+1,:,g);
    for x=2:200
        newm1(x,:,g)=Dm*m(x-1,:,g)+(1-2*Dm)*m(x,:,g)+Dm*m(x+1,:,g);
    end
    newm1(201,:,g)=Dm*m(201-1,:,g)+(1-Dm)*m(201,:,g);
    newm2(:,1,g)=(1-Dm)*newm1(:,1,g)+Dm*newm1(:,1+1,g);
    for y=2:200
        newm2(:,y,g)=Dm*newm1(:,y-1,g)+(1-2*Dm)*newm1(:,y,g)+Dm*newm1(:,y+1,g);
    end
    newm2(:,201,g)=Dm*newm1(:,201-1,g)+(1-Dm)*newm1(:,201,g);
    for x=1:201
        for y=1:201
            newm1(x,y,g)=newm2(x,y,g)+(theta*Cm(x,y,g)-lambda*m(x,y,g))*deltat;
        end
    end
end
%2. ECM ALL SITES
for x=1:201
    for y=1:201
        for g=1:G
            neww(x,y,g)=w(x,y,g)-(alpha1*Cm(x,y,g)+alpha2*m(x,y,g))*w(x,y,g)*deltat;
        end
    end
end
%3. Update microenvironment
for x=1:201
    for y=1:201
        for g=1:G
            if newm1(x,y,g)>1
                m(x,y,g)=1;
            else
                m(x,y,g)=newm1(x,y,g);
            end
        end
    end
end
w=neww;
%4. Update tumor cells in All Sites
%4.1 mesenchymal-like cells
chosen=zeros(201,201,G);
for x=2:200
    for y=2:200
            for g=1:G
                if dead(x,y,g)>0
                    deadage(x,y,g)=deadage(x,y,g)+1;
                    if deadage(x,y,g)>=Td
                           dead(x,y,g)=deadage(x,y,g)-1;
                           deadage(x,y,g)=0;
                    end
                end
                if Cm(x,y,g)>0 && chosen(x,y,g)==0
                    chosen(x,y,g)=1;
                    if agecm(x,y,g)>=Tdm
                        Cm(x,y,g)=Cm(x,y,g)-1;
                        if Cm(x,y,g)==0
                        agecm(x,y,g)=0;
                        else
                        agecm(x,y,g)=1;
                        end
                        dead(x,y,g)=dead(x,y,g)+1;
                        deadage(x,y,g)=deadage(x,y,g)+1;
                        else
                           agecm(x,y,g)=agecm(x,y,g)+1;
                        if agecm(x,y,g)>=Tm
                            if n(x+1,y,g)<4 || n(x-1,y,g)<4 || n(x,y+1,g)<4 || n(x,y-1,g)<4
                                key=floor(rand(1)*4);
                                switch key
                                    case 0
                                        if n(x+1,y,g)<4
                                            Cm(x+1,y,g)=Cm(x+1,y,g)+1;
                                            if dead(x+1,y,g)>0
                                            dead(x+1,y,g)=dead(x+1,y,g)-1;
                                            else
                                            dead(x+1,y,g)=0;
                                            end
                                            if agecm(x+1,y,g)>0
                                            agecm(x+1,y,g)=agecm(x+1,y,g)+1;
                                            else
                                            agecm(x+1,y,g)=1;
                                            end
                                            chosen(x+1,y,g)=1;
                                            key=-1;
                                        end
                                    case 1
                                        if n(x-1,y,g)<4
                                            Cm(x-1,y,g)=Cm(x-1,y,g)+1;
                                            if dead(x-1,y,g)>0
                                            dead(x-1,y,g)=dead(x-1,y,g)-1;
                                            else
                                            dead(x-1,y,g)=0;
                                            end
                                            if agecm(x-1,y,g)>0
                                            agecm(x-1,y,g)=agecm(x-1,y,g)+1;
                                            else
                                            agecm(x-1,y,g)=1;
                                            end
                                            chosen(x-1,y,g)=1;
                                            key=-1;                                    
                                        end
                                    case 2
                                        if n(x,y+1,g)<4
                                            Cm(x,y+1,g)=Cm(x,y+1,g)+1;
                                            if dead(x,y+1,g)>0
                                            dead(x,y+1,g)=dead(x,y+1,g)-1;
                                            else
                                            dead(x,y+1,g)=0;
                                            end
                                            if agecm(x,y+1,g)>0
                                            agecm(x,y+1,g)=agecm(x,y+1,g)+1;
                                            else
                                            agecm(x,y+1,g)=1;
                                            end
                                            chosen(x,y+1,g)=1;
                                            key=-1;                                    
                                        end
                                    case 3
                                        if n(x,y-1,g)<4
                                            Cm(x,y-1,g)=Cm(x,y-1,g)+1;
                                            if dead(x,y-1,g)>0
                                            dead(x,y-1,g)=dead(x,y-1,g)-1;
                                            else
                                            dead(x,y-1,g)=0;
                                            end
                                            if agecm(x,y-1,g)>0
                                            agecm(x,y-1,g)=agecm(x,y-1,g)+1;
                                            else
                                            agecm(x,y-1,g)=1;
                                            end
                                            chosen(x,y-1,g)=1;
                                            key=-1;                                    
                                        end
                                end
                                if key==-1
                                   key=rand(1);
                                   if key>0.1
                                      Cm(x,y,g)=Cm(x,y,g)-1;
                                   end
                                end                        
                             end
                         end
                    end
                 end
            end
    end
end
chosen=zeros(201,201,G);
for x=2:200
    for y=2:200
            for g=1:G
                if Cm(x,y,g)>0 && chosen(x,y,g)==0 
                    P0=1-4*deltat*Dcm/deltax/deltax-deltat*phicm/deltax/deltax*(w(x+1,y,g)+w(x-1,y,g)+w(x,y+1,g)+w(x,y-1,g)-4*w(x,y,g));
                    P1=deltat*Dcm/deltax/deltax-deltat*phicm/4/deltax/deltax*(w(x+1,y,g)-w(x-1,y,g));
                    P2=deltat*Dcm/deltax/deltax+deltat*phicm/4/deltax/deltax*(w(x+1,y,g)-w(x-1,y,g));
                    P3=deltat*Dcm/deltax/deltax-deltat*phicm/4/deltax/deltax*(w(x,y+1,g)-w(x,y-1,g));
                    P4=deltat*Dcm/deltax/deltax+deltat*phicm/4/deltax/deltax*(w(x,y+1,g)-w(x,y-1,g));
                    if P0<0
                        P0=0;
                    end
                    if P1<0
                        P1=0;
                    end
                    if P2<0
                        P2=0;
                    end
                    if P3<0
                        P3=0;
                    end
                    if P4<0
                        P4=0;
                    end
                    PP0=P0;
                    PP1=P0+P1;
                    PP2=P0+P1+P2;
                    PP3=P0+P1+P2+P3;
                    PP4=P0+P1+P2+P3+P4;
                    key=rand(1)*PP4;
                    if key<PP0    
                        chosen(x,y,g)=1;
                    elseif key<PP1
                        if n(x-1,y,g)<4
                        Cm(x-1,y,g)=Cm(x-1,y,g)+1;
                        dead(x-1,y,g)=0;
                        agecm(x-1,y,g)=agecm(x,y,g);
                        chosen(x-1,y,g)=1;
                        Cm(x,y,g)=Cm(x,y,g)-1;
                        end
                    elseif key<PP2
                        if n(x+1,y,g)<4
                        Cm(x+1,y,g)=Cm(x+1,y,g)+1;
                        dead(x+1,y,g)=0;
                        agecm(x+1,y,g)=agecm(x,y,g);
                        chosen(x+1,y,g)=1;
                        Cm(x,y,g)=Cm(x,y,g)-1;
                        end
                    elseif key<PP3
                        if n(x,y-1,g)<4
                        Cm(x,y-1,g)=Cm(x,y-1,g)+1;
                        dead(x,y-1,g)=0;
                        agecm(x,y-1,g)=agecm(x,y,g);
                        chosen(x,y-1,g)=1;
                        Cm(x,y,g)=Cm(x,y,g)-1;    
                        end
                    else
                        if n(x,y+1,g)<4
                        Cm(x,y+1,g)=Cm(x,y+1,g)+1;
                        dead(x,y+1,g)=0;
                        agecm(x,y+1,g)=agecm(x,y,g);
                        chosen(x,y+1,g)=1;
                        Cm(x,y,g)=Cm(x,y,g)-1; 
                        end
                    end
                end
            end
    end
end
%4.1 epithelial-like cells
chosen=zeros(201,201,G);
for x=2:200
    for y=2:200
            for g=1:G
                if Ce(x,y,g)>0 && chosen(x,y,g)==0
                    chosen(x,y,g)=1;
                    if agece(x,y,g)>=Tde
                        Ce(x,y,g)=Ce(x,y,g)-1;
                        if Ce(x,y,g)==0
                        agece(x,y,g)=0;
                        else
                        agece(x,y,g)=1;
                        end
                        dead(x,y,g)=dead(x,y,g)+1;
                        deadage(x,y,g)=deadage(x,y,g)+1;
                        else
                           agece(x,y,g)=agece(x,y,g)+1;
                        if agece(x,y,g)>=Te
                            if n(x+1,y,g)<4 || n(x-1,y,g)<4 || n(x,y+1,g)<4 || n(x,y-1,g)<4
                                key=floor(rand(1)*4);
                                switch key
                                    case 0
                                        if n(x+1,y,g)<4
                                            Ce(x+1,y,g)=Ce(x+1,y,g)+1;
                                            if dead(x+1,y,g)>0
                                            dead(x+1,y,g)=dead(x+1,y,g)-1;
                                            else
                                            dead(x+1,y,g)=0;
                                            end
                                            if agece(x+1,y,g)>0
                                            agece(x+1,y,g)=agece(x+1,y,g)+1;
                                            else
                                            agece(x+1,y,g)=1;
                                            end
                                            chosen(x+1,y,g)=1;
                                            key=-1;
                                        end
                                    case 1
                                        if n(x-1,y,g)<4
                                            Ce(x-1,y,g)=Ce(x-1,y,g)+1;
                                            if dead(x-1,y,g)>0
                                            dead(x-1,y,g)=dead(x-1,y,g)-1;
                                            else
                                            dead(x-1,y,g)=0;
                                            end
                                            if agece(x-1,y,g)>0
                                            agece(x-1,y,g)=agece(x-1,y,g)+1;
                                            else
                                            agece(x-1,y,g)=1;
                                            end
                                            chosen(x-1,y,g)=1;
                                            key=-1;                                    
                                        end
                                    case 2
                                        if n(x,y+1,g)<4
                                            Ce(x,y+1,g)=Ce(x,y+1,g)+1;
                                            if dead(x,y+1,g)>0
                                            dead(x,y+1,g)=dead(x,y+1,g)-1;
                                            else
                                            dead(x,y+1,g)=0;
                                            end
                                            if agece(x,y+1,g)>0
                                            agece(x,y+1,g)=agece(x,y+1,g)+1;
                                            else
                                            agece(x,y+1,g)=1;
                                            end
                                            chosen(x,y+1,g)=1;
                                            key=-1;                                    
                                        end
                                    case 3
                                        if n(x,y-1,g)<4
                                            Ce(x,y-1,g)=Ce(x,y-1,g)+1;
                                            if dead(x,y-1,g)>0
                                            dead(x,y-1,g)=dead(x,y-1,g)-1;
                                            else
                                            dead(x,y-1,g)=0;
                                            end
                                            if agece(x,y-1,g)>0
                                            agece(x,y-1,g)=agece(x,y-1,g)+1;
                                            else
                                            agece(x,y-1,g)=1;
                                            end
                                            chosen(x,y-1,g)=1;
                                            key=-1;                                    
                                        end
                                end
                                if key==-1
                                   key=rand(1);
                                   if key>0.1
                                      Ce(x,y,g)=Ce(x,y,g)-1;
                                   end
                                end                        
                             end
                         end
                    end
                 end
            end
    end
end
chosen=zeros(201,201,G);
for x=2:200
    for y=2:200
            for g=1:G
                if Ce(x,y,g)>0 && chosen(x,y,g)==0 
                    P0=1-4*deltat*Dce/deltax/deltax-deltat*phice/deltax/deltax*(w(x+1,y,g)+w(x-1,y,g)+w(x,y+1,g)+w(x,y-1,g)-4*w(x,y,g));
                    P1=deltat*Dce/deltax/deltax-deltat*phice/4/deltax/deltax*(w(x+1,y,g)-w(x-1,y,g));
                    P2=deltat*Dce/deltax/deltax+deltat*phice/4/deltax/deltax*(w(x+1,y,g)-w(x-1,y,g));
                    P3=deltat*Dce/deltax/deltax-deltat*phice/4/deltax/deltax*(w(x,y+1,g)-w(x,y-1,g));
                    P4=deltat*Dce/deltax/deltax+deltat*phice/4/deltax/deltax*(w(x,y+1,g)-w(x,y-1,g));
                    if P0<0
                        P0=0;
                    end
                    if P1<0
                        P1=0;
                    end
                    if P2<0
                        P2=0;
                    end
                    if P3<0
                        P3=0;
                    end
                    if P4<0
                        P4=0;
                    end
                    PP0=P0;
                    PP1=P0+P1;
                    PP2=P0+P1+P2;
                    PP3=P0+P1+P2+P3;
                    PP4=P0+P1+P2+P3+P4;
                    key=rand(1)*PP4;
                    if key<PP0    
                        chosen(x,y,g)=1;
                    elseif key<PP1
                        if n(x-1,y,g)<4
                        Ce(x-1,y,g)=Ce(x-1,y,g)+1;
                        dead(x-1,y,g)=0;
                        agece(x-1,y,g)=agece(x,y,g);
                        chosen(x-1,y,g)=1;
                        Ce(x,y,g)=Ce(x,y,g)-1;
                        end
                    elseif key<PP2
                        if n(x+1,y,g)<4
                        Ce(x+1,y,g)=Ce(x+1,y,g)+1;
                        dead(x+1,y,g)=0;
                        agece(x+1,y,g)=agece(x,y,g);
                        chosen(x+1,y,g)=1;
                        Ce(x,y,g)=Ce(x,y,g)-1;
                        end
                    elseif key<PP3
                        if n(x,y-1,g)<4
                        Ce(x,y-1,g)=Ce(x,y-1,g)+1;
                        dead(x,y-1,g)=0;
                        agece(x,y-1,g)=agece(x,y,g);
                        chosen(x,y-1,g)=1;
                        Ce(x,y,g)=Ce(x,y,g)-1;    
                        end
                    else
                        if n(x,y+1,g)<4
                        Ce(x,y+1,g)=Ce(x,y+1,g)+1;
                        dead(x,y+1,g)=0;
                        agece(x,y+1,g)=agece(x,y,g);
                        chosen(x,y+1,g)=1;
                        Ce(x,y,g)=Ce(x,y,g)-1; 
                        end
                    end
                end

            end
    end
end
%5. intravasation
for x=1:201
    for y=1:201
            if Cm(x,y,1)>0 && u(x,y,1)==1
                Vcm(x,y)=Vcm(x,y)+Cm(x,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x+1,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x-1,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x,y+1,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x,y-1,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x+1,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x-1,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y+1,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y-1,1);
                Cm(x,y,1)=0;
                Cm(x+1,y,1)=0;
                Cm(x-1,y,1)=0;
                Cm(x,y+1,1)=0;
                Cm(x,y-1,1)=0;
                agecm(x,y,1)=0;
                agecm(x+1,y,1)=0;
                agecm(x-1,y,1)=0;
                agecm(x,y+1,1)=0;
                agecm(x,y-1,1)=0;
                Ce(x,y,1)=0;
                Ce(x+1,y,1)=0;
                Ce(x-1,y,1)=0;
                Ce(x,y+1,1)=0;
                Ce(x,y-1,1)=0;
                agece(x,y,1)=0;
                agece(x+1,y,1)=0;
                agece(x-1,y,1)=0;
                agece(x,y+1,1)=0;
                agece(x,y-1,1)=0;
            end
            if n(x,y,1)>0 && v(x,y,1)==1
                Vcm(x,y)=Vcm(x,y)+Cm(x,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x+1,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x-1,y,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x,y+1,1);
                Vcm(x,y)=Vcm(x,y)+Cm(x,y-1,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x+1,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x-1,y,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y+1,1);
                Vce(x,y)=Vce(x,y)+Ce(x,y-1,1);
                Cm(x,y,1)=0;
                Cm(x+1,y,1)=0;
                Cm(x-1,y,1)=0;
                Cm(x,y+1,1)=0;
                Cm(x,y-1,1)=0;
                agecm(x,y,1)=0;
                agecm(x+1,y,1)=0;
                agecm(x-1,y,1)=0;
                agecm(x,y+1,1)=0;
                agecm(x,y-1,1)=0;
                Ce(x,y,1)=0;
                Ce(x+1,y,1)=0;
                Ce(x-1,y,1)=0;
                Ce(x,y+1,1)=0;
                Ce(x,y-1,1)=0;
                agece(x,y,1)=0;
                agece(x+1,y,1)=0;
                agece(x-1,y,1)=0;
                agece(x,y+1,1)=0;
                agece(x,y-1,1)=0;
            end
    end
end

%6. travel through vessel
for x=1:201
    for y=1:201
        if Vce(x,y)==0 && Vcm(x,y)==0
            agev(x,y)=0;
        else
            agev(x,y)=agev(x,y)+1;
        end
        if Vce(x,y)+Vcm(x,y)>1 && agev(x,y)==floor(Tv/2)
            key=rand(1);
            if key<0.2
                i=randi([1,201]);
                j=randi([1,201]);
                Vcm(i,j)=Vcm(i,j)+1;
                Vcm(x,y)=Vcm(x,y)-1;
            end
        end
        if (Vce(x,y)>0 || Vcm(x,y)>0) && agev(x,y)>=Tv
            key=rand(1);
            if Vce(x,y)+Vcm(x,y)>1
                if key>Pc
                    Vce(x,y)=0;
                    Vcm(x,y)=0;
                    agev(x,y)=0;
                end
            else
                if key>Ps
                    Vce(x,y)=0;
                    Vcm(x,y)=0;
                    agev(x,y)=0;
                end
            end
        end
    end
end
%7. Extravasation 
for x=1:201
    for y=1:201
        if Vce(x,y)+Vcm(x,y)>0 && agev(x,y)>=Tv
            pp1=ebslon3;
            pp2=ebslon3+ebslon2;
            pp3=ebslon3+ebslon2+ebslon1;
            key=rand(1)*pp3;
            if key<pp1
                [i,j]=find(u(:,:,2)==1);
                idx=randi([1,10]);
                i1=i(idx);
                j1=j(idx);
                for num=1:Vce(x,y)+Vcm(x,y)
                    if n(i1,j1,2)<4 && Vcm(x,y)>0
                        Cm(i1,j1,2)=Cm(i1,j1,2)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1+1,j1,2)<4 && Vcm(x,y)>0
                        Cm(i1+1,j1,2)=Cm(i1+1,j1,2)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1-1,j1,2)<4 && Vcm(x,y)>0
                        Cm(i1-1,j1,2)=Cm(i1-1,j1,2)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1+1,2)<4 && Vcm(x,y)>0
                        Cm(i1,j1+1,2)=Cm(i1,j1+1,2)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1-1,2)<4 && Vcm(x,y)>0
                        Cm(i1,j1-1,2)=Cm(i1,j1-1,2)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    end
                    if n(i1,j1,2)<4 && Vce(x,y)>0
                        Ce(i1,j1,2)=Ce(i1,j1,2)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1+1,j1,2)<4 && Vce(x,y)>0
                        Ce(i1+1,j1,2)=Ce(i1+1,j1,2)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1-1,j1,2)<4 && Vce(x,y)>0
                        Ce(i1-1,j1,2)=Ce(i1-1,j1,2)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1+1,2)<4 && Vce(x,y)>0
                        Ce(i1,j1+1,2)=Ce(i1,j1+1,2)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1-1,2)<4 && Vce(x,y)>0
                        Ce(i1,j1-1,2)=Ce(i1,j1-1,2)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    end
                    n(:,:,2)=Cm(:,:,2)+Ce(:,:,2);
                end
                Vce(x,y)=0;
                Vcm(x,y)=0;
                agev(x,y)=0;
            end  
            if key<pp2
                [i,j]=find(u(:,:,3)==1);
                idx=randi([1,10]);
                i1=i(idx);
                j1=j(idx);
                for num=1:Vce(x,y)+Vcm(x,y)
                    if n(i1,j1,3)<4 && Vcm(x,y)>0
                        Cm(i1,j1,3)=Cm(i1,j1,3)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1+1,j1,3)<4 && Vcm(x,y)>0
                        Cm(i1+1,j1,3)=Cm(i1+1,j1,3)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1-1,j1,3)<4 && Vcm(x,y)>0
                        Cm(i1-1,j1,3)=Cm(i1-1,j1,3)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1+1,3)<4 && Vcm(x,y)>0
                        Cm(i1,j1+1,3)=Cm(i1,j1+1,3)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1-1,3)<4 && Vcm(x,y)>0
                        Cm(i1,j1-1,3)=Cm(i1,j1-1,3)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    end
                    if n(i1,j1,3)<4 && Vce(x,y)>0
                        Ce(i1,j1,3)=Ce(i1,j1,3)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1+1,j1,3)<4 && Vce(x,y)>0
                        Ce(i1+1,j1,3)=Ce(i1+1,j1,3)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1-1,j1,3)<4 && Vce(x,y)>0
                        Ce(i1-1,j1,3)=Ce(i1-1,j1,3)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1+1,3)<4 && Vce(x,y)>0
                        Ce(i1,j1+1,3)=Ce(i1,j1+1,3)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1-1,3)<4 && Vce(x,y)>0
                        Ce(i1,j1-1,3)=Ce(i1,j1-1,3)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    end
                    n(:,:,3)=Cm(:,:,3)+Ce(:,:,3);
                end
                Vce(x,y)=0;
                Vcm(x,y)=0;
                agev(x,y)=0;
            end
            if key<pp3
                [i,j]=find(u(:,:,4)==1);
                idx=randi([1,10]);
                i1=i(idx);
                j1=j(idx);
                for num=1:Vce(x,y)+Vcm(x,y)
                    if n(i1,j1,4)<4 && Vcm(x,y)>0
                        Cm(i1,j1,4)=Cm(i1,j1,4)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1+1,j1,4)<4 && Vcm(x,y)>0
                        Cm(i1+1,j1,4)=Cm(i1+1,j1,4)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1-1,j1,4)<4 && Vcm(x,y)>0
                        Cm(i1-1,j1,4)=Cm(i1-1,j1,4)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1+1,4)<4 && Vcm(x,y)>0
                        Cm(i1,j1+1,4)=Cm(i1,j1+1,4)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    elseif n(i1,j1-1,4)<4 && Vcm(x,y)>0
                        Cm(i1,j1-1,4)=Cm(i1,j1-1,4)+1;
                        Vcm(x,y)=Vcm(x,y)-1;
                    end
                    if n(i1,j1,4)<4 && Vce(x,y)>0
                        Ce(i1,j1,4)=Ce(i1,j1,4)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1+1,j1,4)<4 && Vce(x,y)>0
                        Ce(i1+1,j1,4)=Ce(i1+1,j1,4)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1-1,j1,4)<4 && Vce(x,y)>0
                        Ce(i1-1,j1,4)=Ce(i1-1,j1,4)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1+1,4)<4 && Vce(x,y)>0
                        Ce(i1,j1+1,4)=Ce(i1,j1+1,4)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    elseif n(i1,j1-1,4)<4 && Vce(x,y)>0
                        Ce(i1,j1-1,4)=Ce(i1,j1-1,4)+1;
                        Vce(x,y)=Vce(x,y)-1;
                    end
                    n(:,:,4)=Cm(:,:,4)+Ce(:,:,4);
                end
                Vce(x,y)=0;
                Vcm(x,y)=0;
                agev(x,y)=0;
            end             
        end
    end
end
for g=1:G
    n(:,:,g)=Cm(:,:,g)+Ce(:,:,g);
end
end
