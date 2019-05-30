
%Entrees

beta=0.0349; %angle de calage en rad
Vinf=0:0.7:15; %variation de la vitesse du vent en m/s
R=5; %distance min de du centre de rotor en m
r=1.25; %distance min de du centre de rotor en m
B=2; %nombre de pale
om=7.5; %vitesse de rotation en rad/s
eps=10^(-6); %precision
cl=0.6; %coef de portance
cd=0.07; %coef de trainee
rho=1.2; %masse volumique de l'air

%discretisation de la pale en 400elements
ri=linspace(r,R,200); 
h=(5-1.25)/400;


n=length(Vinf);
m=length(ri);

sigma=zeros(n);
Phi=zeros(m,n);
c=zeros(m); %corde
s=zeros(m); %surface

alpha=zeros(m,n); %angle d'attaque
w=zeros(m,n); %vitesse relative
cpe=zeros(m,n); %coef de puissance pour un element
cptot=zeros(1,n); %coef de puissance globale
Pmax=zeros(1,n);
Ptot=zeros(1,n); %puissance totale

%Initialisation
a1=0.3*ones(m,n);
a2=zeros(m,n);
temp=[]; %var temporaire


%calcul élémentaire
for j=1:n
    a1=0.3*ones(m,n);
    a2=zeros(m,n);
    for i=1:m
        c(i)=0.864-0.1016*ri(i);
        s(i)=0.0375*c(i);
        sigma(i)=(B*c(i))/(2*pi*ri(i));
        temp(1)=0;
        temp(2)=0.3;
        while (abs(a1(i,j)-temp(1))>eps && abs(a2(i,j)-temp(2))>eps) 
         Phi(i,j)=atan((Vinf(j)*(1-a1(i,j)))/((om*ri(i))*(1+a2(i,j))));
         w(i,j)=sqrt(((Vinf(j)*(1-a1(i,j)))^2)+(ri(i)*om*(1+a2(i,j)))^2);
         alpha(i,j)=(Phi(i,j)-beta)*57.2958;
         temp(1)=a1(i,j);
         temp(2)=a2(i,j);
        a1(i,j)=(1+((4*(sin(Phi(i,j))^2)/(sigma(i)*(cl*cos(Phi(i,j))+cd*sin(Phi(i,j)))))))^(-1);
        a2(i,j)=(-1+((4*(sin(Phi(i,j))*cos(Phi(i,j)))/(sigma(i)*(-cd*cos(Phi(i,j))+cl*sin(Phi(i,j)))))))^(-1);
        end

        cpe(i,j)=8*((Vinf(j)*R)^(-2))*(om^2)*(1-a1(i,j))*a2(i,j)*(ri(i)^3);
    end
end

ri=transpose(ri);

%calcul de l'integrale des éléments par la méthode de trapezes
for j=1:n
   cptot(1,j)=B*trapz(ri(:,1),cpe(:,j));
end

%cptot(10)=0.12;

for j=1:n
Pmax(j)=2*0.5*pi*rho*(R^2)*(Vinf(j))^2;
end 

%calcul de la puissance totale
for j=1:n
    Ptot(j)=cptot(j)*Pmax(j);
end

 %plot(Vinf,cptot)
 i=0;
 for V=15:0.1:20
     i=i+1;
     Ptot(i+22)=Ptot(22);
     Vinf(i+22)=V;
 end
 

plot(Vinf,Ptot)
xlabel('Vitesse du vent (m/s)')
ylabel('Puissance extraite')