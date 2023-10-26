%% PROYECTO: TESIS DE MAESTRIA
% AUTOR: Juan Pantoja Moyano.
% DETALLE: Creación de archivo de ingreso de datos Perform 3D
% MODELO: Muelle PEND 2:1 - 323M

%% PREAMBULO

clc
clear
close all
%% INGRESO DE VARIABLES

Ly=[1.5,2.3,2.3,3.05,3.05,3.05,3.05,3.05,3.05,2.3,2.3,1.5,1.5,1.1,0.4];     %Vector de espaciamientos en el sentido Y
Lx=2.5;                                                                     %Separación tipica en x
Lz=1.75;                                                                    %Altura donde esta la placa.
Lxi=[-4,-2.5,2.5,4];                                                        %Vector de espaciamiento Inicial
N=43;                                                                       %Numero de porticos

%% DATOS DE SECCIONES PILOTES, PLACA Y TERRENO
NvTERR=-0.2457;         %[m] Altura en Z al nivel de terreno en Y=0
Hlosa=0.5;              %[m] espesor de la placa
LongDes=2;              %[m] Longitud de desarrollo
Pend=2;                 %Pendiente del talud %:1

% Dimensiones Viga Cabezal y grua

Dim_vc=[1.5 0.7];   %[b h] en m
Dim_vg=[1.5 1.0];   %[b h] en m
rccto=23.536;       %[kN/m^3]

% Control de consola

ctr1=1;                 %(0=EJE 1, 1=TODOS EJES) Se crean todos los ejes de pilotes o solo el primer eje?
ctr2=1;                 %(0=ARCHIVOS ADICIONALES DE CARGA VERTICAL, 1=FRAMES + RESUMEN MASA Y CARGA VERTICAL) OPCION 1: Se crean archivos  de elementos ademas de un archivo resumido de Masa y cargas verticales;  OPCION 2: Se crean archivos individuales de carga vertical.

% DATOS DE SEPARACIÓN DE ELEMENTOS LINK
Link_type ={'K126-1','K126-2','K126-3','K246-1','K390-1'};
LinkLim=[126/3,126/3,126/3,246,390]*0.0254;        %Intervalos entre cambios de tipo de link
LinkSpac=[12,12,12,24,48]*0.0254;                   %Espaciamiento entre elementos Link para cada intervalo

% DATOS DE SEPARACIÓN DE ELEMENTOS PILOTE

PileDiv=[4,4,4,4,4];                   %Numero de elementos Pile entre cada tipo de link


%% CARGAS APLICAR
% D_PLACA
D_PLACA=rccto*0.35; %Placa de 35cm
% D_CARPETA
D_CARPETA=3;        % kN/m^2
% D_SISTEMA_BERTHING
D_SISBER=52;         % kN
Nod_SISBER=3:6:N*3;
% D_EQUIPO
D_EQUIP=1.50;       % kN/m^2
% L_UNIFORME_C
L_UNIFC=[38,19]*0.1;    % kN/m^2
Fil_UNIC=[1,1,1,1,1,1,1,1,1,1,1,1,2,2,2];     %Vector de tipo de carga en el sentido Y
Fil_Divi=[13];                                %Vector de filas que tienen area aferente dividida en dos cargas  
% MASA DE GRUA (Distribuida sobre Viga Grua)
D_CRAINE=2;         %kN/m

%% MASAS DE PILOTES A APLICAR

MasaPilvc=[2.246,3.124,4.002,4.88]*9.81;              %Masa en los nodos desde A,B,....F
MasaPilvg=[1.368,5.758]*9.81;
Col_Masvc=3:3:N*3;                              %Columnas en la matriz de coordenadas donde estan los porticos sentido Viga Cabezal
Col_Masvg=2:N*3+1;                              %Columnas en la matriz de coordenadas donde estan los porticos sentido Viga Grua
Fil_Masvc=[4,6,8,10]; 
Fil_Masvg=[1,13];
Apoyos_Pilvc=[-17,-19,-21,-23];                 %Cota de apoyos para los pilotes en el eje de vigas cabezales
Apoyos_Pilvg=[-15,-26];                         %Cota de apoyos para los pilotes en el eje de vigas grua
OrdenPile=[2,1,1,1,1,2];                        %Se coloca en orden desde A,B... que tipo de pilote es (VC o VG)

%% GENERADOR DE ARCHIVOS LINKS.

for i=1:length(MasaPilvc)
    
    Link_Vc(1,i)=NvTERR-sum(Ly(1:Fil_Masvc(i)-1))/Pend-0.5*LinkSpac(1);
    H=Link_Vc(1,i);
    cont=2;
    Limit=1;
    contLk=1;
    
    while H>Apoyos_Pilvc(i)
        H=H-LinkSpac(Limit);
        if abs(H-Link_Vc(1,i)-0.5*LinkSpac(1))-1e-5<=sum(LinkLim(1:Limit))
            Link_Vc(cont,i)=H;
            cont=cont+1;
            contLk=contLk+1;
        elseif Limit<length(LinkLim)
            H=Link_Vc(cont-1,i)-0.5*(LinkSpac(Limit)+LinkSpac(Limit+1));
            Link_Vc(cont,i)=H;
            NlinkVc(Limit,i)=contLk;
            Limit=Limit+1;
            cont=cont+1;
            contLk=1;
        else
            H=Link_Vc(cont-1,i)-LinkSpac(Limit);
            Link_Vc(cont,i)=H;
            cont=cont+1;
            contLk=contLk+1;
        end
    end
    Link_Vc(cont-1,i)=Apoyos_Pilvc(i);
    NlinkVc(Limit,i)=contLk-1;
    contLk=1;
end

for i=1:length(MasaPilvg)
    
    Link_Vg(1,i)=NvTERR-sum(Ly(1:Fil_Masvg(i)-1))/Pend-0.5*LinkSpac(1);
    H=Link_Vg(1,i);
    cont=2;
    Limit=1;
    contLk=1;
    
    while H>Apoyos_Pilvg(i)
        H=H-LinkSpac(Limit);
        if abs(H-Link_Vg(1,i)-0.5*LinkSpac(1))-1e-5<=sum(LinkLim(1:Limit))
            Link_Vg(cont,i)=H;
            cont=cont+1;
            contLk=contLk+1;
        elseif Limit<length(LinkLim)
            H=Link_Vg(cont-1,i)-0.5*(LinkSpac(Limit)+LinkSpac(Limit+1));
            Link_Vg(cont,i)=H;
            NlinkVg(Limit,i)=contLk;
            Limit=Limit+1;
            cont=cont+1;
            contLk=1;
        else
            H=Link_Vg(cont-1,i)-LinkSpac(Limit);
            Link_Vg(cont,i)=H;
            cont=cont+1;
            contLk=contLk+1;
        end
    end
    Link_Vg(cont-1,i)=Apoyos_Pilvg(i);
    NlinkVg(Limit,i)=contLk-1;
    contLk=1;
end

%% GENERADOR DE ARCHIVOS PILOTES..
    P=3;
for i=1:length(MasaPilvc)
    
    Pile_Vc(1,i)=Lz;
    Pile_Vc(2,i)=Lz-LongDes;
    cont=P;
    Limit=1;
    
    for k=1:sum(NlinkVc(:,i))-1
        if Link_Vc(k+1,i)~=0
            if k<= sum(NlinkVc(1:Limit,i))-1
            elseif Limit<length(LinkLim)
                    Limit=Limit+1;
            end
            Space=linspace(Link_Vc(k,i),Link_Vc(k+1,i),PileDiv(Limit)+1);
            for l=1:length(Space)-1
                Pile_Vc(cont,i)=Space(l);
                cont=cont+1;
            end
        end
    end
    Pile_Vc(cont,i)=Space(end);
    Pile_Vc(cont+1,i)=Apoyos_Pilvc(i);
end

for i=1:length(MasaPilvg)
    
    Pile_Vg(1,i)=Lz;
    Pile_Vg(2,i)=Lz-LongDes;
    cont=P;
    Limit=1;
    
    for k=1:sum(NlinkVg(:,i))-1
        if Link_Vg(k+1,i)~=0
            if k<= sum(NlinkVg(1:Limit,i))-1
            elseif Limit<length(LinkLim)
                    Limit=Limit+1;
            end
            Space=linspace(Link_Vg(k,i),Link_Vg(k+1,i),PileDiv(Limit)+1);
            for l=1:length(Space)-1
                Pile_Vg(cont,i)=Space(l);
                cont=cont+1;
            end
        end
    end
    Pile_Vg(cont,i)=Space(end);
    Pile_Vg(cont+1,i)=Apoyos_Pilvg(i);
end

 %% CALCULO CORDENADAS DE PUNTOS DE PLACA

Coord=zeros(length(Ly),N*3+2,3);

for i=1:length(Ly)+1
    for j=0:(N-1)*3+1
        if j==0
            Coord(i,1,1)=Lxi(1);
            Coord(i,2,1)=Lxi(2);
            Coord(i,N*3+1,1)=(N-1)*3*Lx+Lxi(3);
            Coord(i,N*3+2,1)=(N-1)*3*Lx+Lxi(4);
        else
        Coord(i,j+2,1)=(j-1)*Lx;
        end
       
    end
end

Var=0;
for i=1:length(Ly)
    for j=1:N*3+2
        Coord(i+1,j,2)=Ly(i)+Var;
        Coord(i,j,3)=Lz;
    end
    Var=Var+Ly(i);
end

for i=1:length(Ly)+1
    for j=1:N*3+2
        Coord(i,j,3)=Lz;
    end
end

%% CREAR MATRIZ DE PUNTOS POR ELEMENTO
MatPuntos=zeros(N*3*length(Ly),12);

ContFil=0;
for i=1:length(Ly)          %Marca la Fila en coord
    for j=1:N*3+1       %Marca la columna en coord
        ContFil=ContFil+1;
        ContCol=1;
        for k=0:1           %Marca el avance vertical (KL)
            for l=0:1       %Marca el avance horizontal (IJ)
                for m=1:3   %Extra los 3 datos de cordenada
                    
                    MatPuntos(ContFil,ContCol)=Coord(i+k,j+l,m);
                    ContCol=ContCol+1;
                end
            end
        end
    end
end

%% CALCULO DE AREAS AFERENTES PARA CARGAS EN LOS NODOS
AreasAf=zeros(length(Ly)+1,N*3+2);
for i=1:length(Ly)+1
    for j=1:N*3+2
        if j>2
            if j<=N*3
                if i==1
                    AreasAf(i,j)=Ly(i)*Lx/2;
                elseif i==length(Ly)+1
                    AreasAf(i,j)=Ly(i-1)*Lx/2;
                else
                    AreasAf(i,j)=(Ly(i)+Ly(i-1))*Lx/2;
                end
                     
            end
        end
    end
end
   
for i=1:length(Ly)+1
    if i==1
        AreasAf(i,1)=Ly(i)*(abs(Lxi(1))-abs(Lxi(2)))/4;
        AreasAf(i,2)=Ly(i)*(abs(Lxi(1)))/4;
        AreasAf(i,N*3+1)=Ly(i)*(abs(Lxi(4)))/4;
        AreasAf(i,N*3+2)=Ly(i)*(abs(Lxi(4))-abs(Lxi(3)))/4;
    elseif i==length(Ly)+1 
        AreasAf(i,1)=Ly(i-1)*(abs(Lxi(1))-abs(Lxi(2)))/4;
        AreasAf(i,2)=Ly(i-1)*(abs(Lxi(1)))/4;
        AreasAf(i,N*3+1)=Ly(i-1)*(abs(Lxi(4)))/4;
        AreasAf(i,N*3+2)=Ly(i-1)*(abs(Lxi(4))-abs(Lxi(3)))/4;
    else
        AreasAf(i,1)=(Ly(i)+Ly(i-1))*(abs(Lxi(1))-abs(Lxi(2)))/4;
        AreasAf(i,2)=(Ly(i)+Ly(i-1))*(abs(Lxi(1)))/4;
        AreasAf(i,N*3+1)=(Ly(i)+Ly(i-1))*(abs(Lxi(4)))/4;
        AreasAf(i,N*3+2)=(Ly(i)+Ly(i-1))*(abs(Lxi(4))-abs(Lxi(3)))/4;
             
    end
        
end

AreasAfAd=zeros(length(Ly)+1,N*3+2,2);
for i=1:length(Fil_Divi)
    for j=1:N*3+2
        if j>2
            if j<=N*3
                AreasAfAd(Fil_Divi(i),j,1)=Ly(Fil_Divi(i))*Lx/2;
                AreasAfAd(Fil_Divi(i),j,2)=Ly(Fil_Divi(i)-1)*Lx/2;
            end
        end
    end
end

for i=1:length(Fil_Divi)
        AreasAfAd(Fil_Divi(i),1,1)=Ly(Fil_Divi(i))*(abs(Lxi(1))-abs(Lxi(2)))/4;
        AreasAfAd(Fil_Divi(i),1,2)=Ly(Fil_Divi(i)-1)*(abs(Lxi(1))-abs(Lxi(2)))/4;
        AreasAfAd(Fil_Divi(i),2,1)=Ly(Fil_Divi(i))*(abs(Lxi(1)))/4;
        AreasAfAd(Fil_Divi(i),2,2)=Ly(Fil_Divi(i)-1)*(abs(Lxi(1)))/4;
        AreasAfAd(Fil_Divi(i),N*3+1,1)=Ly(Fil_Divi(i))*(abs(Lxi(4)))/4;
        AreasAfAd(Fil_Divi(i),N*3+1,2)=Ly(Fil_Divi(i)-1)*(abs(Lxi(4)))/4;
        AreasAfAd(Fil_Divi(i),N*3+2,1)=Ly(Fil_Divi(i))*(abs(Lxi(4))-abs(Lxi(3)))/4;
        AreasAfAd(Fil_Divi(i),N*3+2,2)=Ly(Fil_Divi(i)-1)*(abs(Lxi(4))-abs(Lxi(3)))/4;
end
 

Areatot=0;
for i=1:length(Ly)+1
    Areatot=Areatot+sum(AreasAf(i,:));
end

%% CALCULO LONGITUDES AFERENTES DE VIGAS CABEZALES Y DE GRUA

% VIGA GRUA
LAfer_Vg=zeros(length(Ly)+1,N*3+2);

for i=1:length(Fil_Masvg)
    for j=3:N*3
        if j>2
            if j<=N*3
                LAfer_Vg(Fil_Masvg(i),j)=Lx;
            end
        end
    end
end
   
for i=1:length(Fil_Masvg)
    LAfer_Vg(Fil_Masvg(i),1)=(abs(Lxi(1))-abs(Lxi(2)))/2;
    LAfer_Vg(Fil_Masvg(i),2)=(abs(Lxi(1)))/2;
    LAfer_Vg(Fil_Masvg(i),N*3+1)=(abs(Lxi(4)))/2;
    LAfer_Vg(Fil_Masvg(i),N*3+2)=(abs(Lxi(4))-abs(Lxi(3)))/2;
end


% VIGA CABEZAL

LAfer_Vc=zeros(length(Ly)+1,N*3+2);

for i=1:length(Ly)+1
    for j=1:length(Col_Masvc)
        if i==1
            LAfer_Vc(i,Col_Masvc(j))=Ly(i)/2;
        elseif i==length(Ly)+1
            LAfer_Vc(i,Col_Masvc(j))=Ly(i-1)/2;
        else
            LAfer_Vc(i,Col_Masvc(j))=(Ly(i)+Ly(i-1))/2;
        end
    end
end

for i=1:length(Fil_Masvg)
    for j=1:length(Col_Masvc)
        LAfer_Vc(Fil_Masvg(i),Col_Masvc(j))=LAfer_Vc(Fil_Masvg(i),Col_Masvc(j))-Dim_vg(1);
        if  LAfer_Vc(Fil_Masvg(i),Col_Masvc(j))<0
            LAfer_Vc(Fil_Masvg(i),Col_Masvc(j))=0;
        end
    end
end

%% CALCULO DE MATRICES DE CARGAS
% D_PLACA
D_MPLACA=D_PLACA*AreasAf;

D_PLACTOT=0;
for i=1:length(Ly)+1
    D_PLACTOT=D_PLACTOT+sum(D_MPLACA(i,:));
end

% D_CARPETA
D_MCARPETA=D_CARPETA*AreasAf;

D_CARPTOT=0;
for i=1:length(Ly)+1
    D_CARPTOT=D_CARPTOT+sum(D_MCARPETA(i,:));
end

% D_SISTEMA_BERTHING

D_MSISBER=zeros(length(Ly)+1,N*3+2);

for i=1:length(Nod_SISBER)
    D_MSISBER(length(Ly)+1,Nod_SISBER(i))=1;
end
    
% D_EQUIPO
D_MEQUIP=D_EQUIP*AreasAf;

D_MEQUIPTOT=0;
for i=1:length(Ly)+1
    D_MEQUIPTOT=D_MEQUIPTOT+sum(D_MEQUIP(i,:));
end

% L_UNIFORME_C

L_MUNIFC=zeros(length(Ly)+1,N*3+2);
for i=1:length(Ly)+1
    if i==length(Ly)+1
        L_MUNIFC(i,:)=AreasAf(i,:)*L_UNIFC(Fil_UNIC(i-1));
    else
        L_MUNIFC(i,:)=AreasAf(i,:)*L_UNIFC(Fil_UNIC(i));
    end
end

for i=1:length(Fil_Divi)
    L_MUNIFC(Fil_Divi(i),:)=AreasAfAd(Fil_Divi(i),:,1)*L_UNIFC(Fil_UNIC(Fil_Divi(i)))+AreasAfAd(Fil_Divi(i),:,2)*L_UNIFC(Fil_UNIC(Fil_Divi(i)-1));
end

L_MUNIFCTOT=0;
for i=1:length(Ly)+1
    L_MUNIFCTOT=L_MUNIFCTOT+sum(L_MUNIFC(i,:));
end

% D_VIGAS
D_MVIGAS=Dim_vc(1)*Dim_vc(2)*rccto*LAfer_Vc+Dim_vg(1)*Dim_vg(2)*rccto*LAfer_Vg;

D_MVIGASTOT=0;
for i=1:length(Ly)+1
    D_MVIGASTOT=D_MVIGASTOT+sum(D_MVIGAS(i,:));
end

% D_CRAINE
D_MCRAINE=D_CRAINE*LAfer_Vg;

D_MCRAINETOT=0;
for i=1:length(Ly)+1
    D_MCRAINETOT=D_MCRAINETOT+sum(D_MCRAINE(i,:));
end
%% CREACIÓN MATRICES DE MASAS

% D_PLACA

M_PLACA=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_PLACA(ContFil,m)=Coord(i,j,m);
        end
        M_PLACA(ContFil,4)=D_MPLACA(i,j);                         %MASA EN H1
        M_PLACA(ContFil,5)=D_MPLACA(i,j);                         %MASA EN H2
        M_PLACA(ContFil,6)=0;                                     %MASA EN H3
        M_PLACA(ContFil,7)=0;                         %MOMENTO EN H2
        M_PLACA(ContFil,8)=0;                         %MOMENTO EN H2
        M_PLACA(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% D_CARPETA

M_CARPETA=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_CARPETA(ContFil,m)=Coord(i,j,m);
        end
        M_CARPETA(ContFil,4)=D_MCARPETA(i,j);                         %MASA EN H1
        M_CARPETA(ContFil,5)=D_MCARPETA(i,j);                         %MASA EN H2
        M_CARPETA(ContFil,6)=0;                                       %MASA EN H3
        M_CARPETA(ContFil,7)=0;                         %MOMENTO EN H2
        M_CARPETA(ContFil,8)=0;                         %MOMENTO EN H2
        M_CARPETA(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% D_SISTEMA BERTHING

M_SISBER=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_SISBER(ContFil,m)=Coord(i,j,m);
        end
        M_SISBER(ContFil,4)=D_SISBER*D_MSISBER(i,j);                         %MASA EN H1
        M_SISBER(ContFil,5)=D_SISBER*D_MSISBER(i,j);                         %MASA EN H2
        M_SISBER(ContFil,6)=0;                                               %MASA EN H3
        M_SISBER(ContFil,7)=0;                         %MASA EN H2
        M_SISBER(ContFil,8)=0;                         %MASA EN H2
        M_SISBER(ContFil,9)=0;                         %MASA EN H3
        
    end
end

% D_EQUIPO

M_EQUIPO=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_EQUIPO(ContFil,m)=Coord(i,j,m);
        end
        M_EQUIPO(ContFil,4)=D_MEQUIP(i,j);                         %MASA EN H1
        M_EQUIPO(ContFil,5)=D_MEQUIP(i,j);                         %MASA EN H2
        M_EQUIPO(ContFil,6)=0;                                     %MASA EN H3
        M_EQUIPO(ContFil,7)=0;                         %MOMENTO EN H2
        M_EQUIPO(ContFil,8)=0;                         %MOMENTO EN H2
        M_EQUIPO(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% L_UNIFORME_C

M_LUNIFC=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_LUNIFC(ContFil,m)=Coord(i,j,m);
        end
        M_LUNIFC(ContFil,4)=L_MUNIFC(i,j);                         %MASA EN H1
        M_LUNIFC(ContFil,5)=L_MUNIFC(i,j);                         %MASA EN H2
        M_LUNIFC(ContFil,6)=0;                                     %MASA EN H3
        M_LUNIFC(ContFil,7)=0;                         %MOMENTO EN H2
        M_LUNIFC(ContFil,8)=0;                         %MOMENTO EN H2
        M_LUNIFC(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% D_VIGAS

M_VIGAS=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_VIGAS(ContFil,m)=Coord(i,j,m);
        end
        M_VIGAS(ContFil,4)=D_MVIGAS(i,j);                         %MASA EN H1
        M_VIGAS(ContFil,5)=D_MVIGAS(i,j);                         %MASA EN H2
        M_VIGAS(ContFil,6)=0;                                     %MASA EN H3
        M_VIGAS(ContFil,7)=0;                         %MOMENTO EN H2
        M_VIGAS(ContFil,8)=0;                         %MOMENTO EN H2
        M_VIGAS(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% D_CRAINE

M_CRAINE=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_CRAINE(ContFil,m)=Coord(i,j,m);
        end
        M_CRAINE(ContFil,4)=D_MCRAINE(i,j);                         %MASA EN H1
        M_CRAINE(ContFil,5)=D_MCRAINE(i,j);                         %MASA EN H2
        M_CRAINE(ContFil,6)=0;                                     %MASA EN H3
        M_CRAINE(ContFil,7)=0;                         %MOMENTO EN H2
        M_CRAINE(ContFil,8)=0;                         %MOMENTO EN H2
        M_CRAINE(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

% MASAS PILOTES

MASA=zeros(length(Ly)+1,N*3+2);

for i=1:length(Col_Masvc)
    for j=1:length(Fil_Masvc)
        MASA(Fil_Masvc(j),Col_Masvc(i))=MasaPilvc(j);
    end
end

for i=1:length(Col_Masvg)
    for j=1:length(Fil_Masvg)
        MASA(Fil_Masvg(j),Col_Masvg(i))=MasaPilvg(j);
    end
end


%% CREACIÓN MATRIZ DE MASAS DE PILOTES

M_MPILO=zeros(N*3*length(Ly),9);

ContFil=0;
for i=1:length(Ly)+1            %Marca la Fila en coord
    for j=1:N*3+2               %Marca la columna en coord
        ContFil=ContFil+1;
        for m=1:3               %Extra los 3 datos de cordenada
            M_MPILO(ContFil,m)=Coord(i,j,m);
        end
        M_MPILO(ContFil,4)=MASA(i,j);                         %MASA EN H1
        M_MPILO(ContFil,5)=MASA(i,j);                         %MASA EN H2
        M_MPILO(ContFil,6)=0;                                     %MASA EN H3
        M_MPILO(ContFil,7)=0;                         %MOMENTO EN H2
        M_MPILO(ContFil,8)=0;                         %MOMENTO EN H2
        M_MPILO(ContFil,9)=0;                         %MOMENTO EN H3
        
    end
end

%% GENERACIÓN ARCHIVO DE FRAMES PILOTES

switch ctr1
    case 0
        m=1;
        n=2;
        o=2;
    case 1
        m=length(Col_Masvc);
        n=1;
        o=length(Col_Masvg);
end

% GENERADOR DE ELEMENTOS PILE - CONEXION

Exp_conex=zeros(length(MasaPilvc)+length(MasaPilvg),1);

ContFil=1;
Cont_Masvc=1;
Cont_Masvg=1;        

for i=1:length(MasaPilvc)+length(MasaPilvg)
    ContF=1;
    switch OrdenPile(i)
        case 1 % El pilote esta en el eje de la viga Viga Cabezal
            for j=1:m                           
                for l=1:1          %Marca Fila en matrix de coordenada Z 
                    if Pile_Vc(l+1,Cont_Masvc)~=0  
                            F_CONEX(ContFil,1)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                            F_CONEX(ContFil,2)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                            F_CONEX(ContFil,3)=Pile_Vc(l+1,Cont_Masvc);
                            F_CONEX(ContFil,4)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                            F_CONEX(ContFil,5)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                            F_CONEX(ContFil,6)=Pile_Vc(l,Cont_Masvc);
                            F_CONEX(ContFil,7)=2;                                               %DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            F_CONEX(ContFil,8)=0;                                               %ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                            ContFil=ContFil+1;
                            ContF=ContF+1;

                    end
                end
            end
            Cont_Masvc=Cont_Masvc+1;

        case 2  % El pilote esta en el eje de la viga Viga de grua
            for j=n:o %Cambio para generar todos los pilotes
                for l=1:1          %Marca Fila en matrix de coordenada Z 
                    if Pile_Vg(l+1,Cont_Masvg)~=0  
                            F_CONEX(ContFil,1)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                            F_CONEX(ContFil,2)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                            F_CONEX(ContFil,3)=Pile_Vg(l+1,Cont_Masvg);
                            F_CONEX(ContFil,4)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                            F_CONEX(ContFil,5)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                            F_CONEX(ContFil,6)=Pile_Vg(l,Cont_Masvg);
                            F_CONEX(ContFil,7)=2;                                               %DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            F_CONEX(ContFil,8)=0;                                               %ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                            ContFil=ContFil+1;
                            ContF=ContF+1;

                    end
                end
            end
            Cont_Masvg=Cont_Masvg+1;
    end
    Exp_conex(i,1)=ContF-1;
end


% GENERADOR DE ELEMENTOS PILE - PILOTE
Exp_pile=zeros(length(MasaPilvc)+length(MasaPilvg),1);

ContFil=1;
Cont_Masvc=1;
Cont_Masvg=1;        

for i=1:length(MasaPilvc)+length(MasaPilvg)
    ContF=1;
    switch OrdenPile(i)
        case 1 % El pilote esta en el eje de la viga Viga Cabezal
            for j=1:m                           
                for l=2:length(Pile_Vc(:,1))-1          %Marca Fila en matrix de coordenada Z 
                    if Pile_Vc(l+1,Cont_Masvc)~=0  
                            F_PILE(ContFil,1)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                            F_PILE(ContFil,2)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                            F_PILE(ContFil,3)=Pile_Vc(l+1,Cont_Masvc);
                            F_PILE(ContFil,4)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                            F_PILE(ContFil,5)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                            F_PILE(ContFil,6)=Pile_Vc(l,Cont_Masvc);
                            F_PILE(ContFil,7)=2;                                               %DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            F_PILE(ContFil,8)=0;                                               %ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                            ContFil=ContFil+1;
                            ContF=ContF+1;

                    end
                end
            end
            Cont_Masvc=Cont_Masvc+1;

        case 2  % El pilote esta en el eje de la viga Viga de grua
            for j=n:o %Cambio para generar todos los pilotes
                for l=2:length(Pile_Vg(:,1))-1          %Marca Fila en matrix de coordenada Z 
                    if Pile_Vg(l+1,Cont_Masvg)~=0  
                            F_PILE(ContFil,1)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                            F_PILE(ContFil,2)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                            F_PILE(ContFil,3)=Pile_Vg(l+1,Cont_Masvg);
                            F_PILE(ContFil,4)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                            F_PILE(ContFil,5)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                            F_PILE(ContFil,6)=Pile_Vg(l,Cont_Masvg);
                            F_PILE(ContFil,7)=2;                                               %DIRECCIÓN DEL EJE LOCAL 2 PARALELO AL $$
                            F_PILE(ContFil,8)=0;                                               %ROTACION DEL EJE LOCAL 2 ALRREDEDOR DEL EJE LOCAL 1
                            ContFil=ContFil+1;
                            ContF=ContF+1;

                    end
                end
            end
            Cont_Masvg=Cont_Masvg+1;
    end
    Exp_pile(i,1)=ContF-1;
end

%% GENERACIÓN ARCHIVO DE FRAMES LINK
Exp_link=zeros(length(LinkLim),1);
ContFil=1;
 
for k=1:length(LinkLim)
    Cont_Masvc=1;
    Cont_Masvg=1;
    ContL=0;
    for i=1:length(MasaPilvc)+length(MasaPilvg)
 
        switch OrdenPile(i)
            case 1 % El pilote esta en el eje de la viga Viga Cabezal
                for j=1:m   %Cambio para generar todos los pilotes                          %Marca la columna en coord
                    if k==1
                        Linf=1;
                    else
                        Linf=sum(NlinkVc(1:k-1,Cont_Masvc))+1;
                    end
                    
                    for l=Linf:sum(NlinkVc(1:k,Cont_Masvc))       %Marca Fila en matrix de coordenada Z 
                        if Link_Vc(l,Cont_Masvc)~=0  
                                F_LINKH(ContFil,1)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1)-1;
                                F_LINKH(ContFil,2)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                                F_LINKH(ContFil,3)=Link_Vc(l,Cont_Masvc);
                                F_LINKH(ContFil,4)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                                F_LINKH(ContFil,5)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                                F_LINKH(ContFil,6)=Link_Vc(l,Cont_Masvc);
                                F_LINKH(ContFil+1,1)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                                F_LINKH(ContFil+1,2)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                                F_LINKH(ContFil+1,3)=Link_Vc(l,Cont_Masvc);
                                F_LINKH(ContFil+1,4)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                                F_LINKH(ContFil+1,5)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2)+1;
                                F_LINKH(ContFil+1,6)=Link_Vc(l,Cont_Masvc);
                                
                                ContFil=ContFil+2;
                                ContL=ContL+2; 
                                
                        end
                    end
                end
                Cont_Masvc=Cont_Masvc+1;
                
                
 
            case 2  % El pilote esta en el eje de la viga Viga de grua
 
                for j=n:o   %Cambio para generar todos los pilotes                          %Marca la columna en coord
                    if k==1
                        Linf=1;
                    else
                        Linf=sum(NlinkVg(1:k-1,Cont_Masvg))+1;
                    end
                    
                    for l=Linf:sum(NlinkVg(1:k,Cont_Masvg))       %Marca Fila en matrix de coordenada Z 
                        if Link_Vg(l,Cont_Masvg)~=0  
                                F_LINKH(ContFil,1)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1)-1;
                                F_LINKH(ContFil,2)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                                F_LINKH(ContFil,3)=Link_Vg(l,Cont_Masvg);
                                F_LINKH(ContFil,4)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                                F_LINKH(ContFil,5)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                                F_LINKH(ContFil,6)=Link_Vg(l,Cont_Masvg);
                                F_LINKH(ContFil+1,1)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                                F_LINKH(ContFil+1,2)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                                F_LINKH(ContFil+1,3)=Link_Vg(l,Cont_Masvg);
                                F_LINKH(ContFil+1,4)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                                F_LINKH(ContFil+1,5)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2)+1;
                                F_LINKH(ContFil+1,6)=Link_Vg(l,Cont_Masvg);
                                
                                ContFil=ContFil+2;
                                ContL=ContL+2; 
                                
                        end
                    end
                end
                Cont_Masvg=Cont_Masvg+1;
        end
    end
    Exp_link(k,1)=ContL;
end


%% GENERACION ARCHIVO DE APOYOS (LINKS Y PILOTES)

ContFil=1;
Cont_Masvc=1;
Cont_Masvg=1;

for i=1:length(MasaPilvc)+length(MasaPilvg)
    
    switch OrdenPile(i)
        case 1 % El pilote esta en el eje de la viga Viga Cabezal
            for j=1:m                           
                    if Pile_Vc(l+1,Cont_Masvc)~=0  
                            N_APOYO(ContFil,1)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),1);
                            N_APOYO(ContFil,2)=Coord(Fil_Masvc(Cont_Masvc),Col_Masvc(j),2);
                            N_APOYO(ContFil,3)=Apoyos_Pilvc(Cont_Masvc);
                            ContFil=ContFil+1;
                    end
            end
            Cont_Masvc=Cont_Masvc+1;

        case 2  % El pilote esta en el eje de la viga Viga de grua
            for j=n:o %Cambio para generar todos los pilotes
                    if Pile_Vg(l+1,Cont_Masvg)~=0  
                            N_APOYO(ContFil,1)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),1);
                            N_APOYO(ContFil,2)=Coord(Fil_Masvg(Cont_Masvg),Col_Masvg(j),2);
                            N_APOYO(ContFil,3)=Apoyos_Pilvg(Cont_Masvg);
                            ContFil=ContFil+1;
                    end
            end
            Cont_Masvg=Cont_Masvg+1;
    end
end

Exp_apoy=[length(N_APOYO(:,1)),length(F_LINKH(:,1))/2,length(F_LINKH(:,1))/2];

N_APOYO=[N_APOYO;F_LINKH(1:2:end,1:3);F_LINKH(2:2:end,4:6)];

%% GENERACIÓN ARCHIVO DE FRAMES VIGAS CABEZALES

ContFil=1;
 
for i=1:length(Col_Masvc)      %Marca Fila en matrix de coordenada Z
    for j=1:length(Coord(:,i))-1
        F_VCAB(ContFil,1)=Coord(j,Col_Masvc(i),1);
        F_VCAB(ContFil,2)=Coord(j,Col_Masvc(i),2);
        F_VCAB(ContFil,3)=Coord(j,Col_Masvc(i),3);
        F_VCAB(ContFil,4)=Coord(j+1,Col_Masvc(i),1);
        F_VCAB(ContFil,5)=Coord(j+1,Col_Masvc(i),2);
        F_VCAB(ContFil,6)=Coord(j+1,Col_Masvc(i),3);
        F_VCAB(ContFil,7)=3;
        F_VCAB(ContFil,8)=0;
        
        
        ContFil=ContFil+1;
    end
end

%% GENERACIÓN ARCHIVO DE FRAMES VIGAS DE GRUA

ContFil=1;
 
for i=1:length(Fil_Masvg)      %Marca Fila en matrix de coordenada Z
    for j=1:length(Coord(i,:,1))-1
        F_VGRUA(ContFil,1)=Coord(Fil_Masvg(i),j,1);
        F_VGRUA(ContFil,2)=Coord(Fil_Masvg(i),j,2);
        F_VGRUA(ContFil,3)=Coord(Fil_Masvg(i),j,3);
        F_VGRUA(ContFil,4)=Coord(Fil_Masvg(i),j+1,1);
        F_VGRUA(ContFil,5)=Coord(Fil_Masvg(i),j+1,2);
        F_VGRUA(ContFil,6)=Coord(Fil_Masvg(i),j+1,3);
        F_VGRUA(ContFil,7)=3;
        F_VGRUA(ContFil,8)=0;
        ContFil=ContFil+1;
    end
end

%% GENERACIÓN ARCHIVO DE FRAMES VIGAS DE BORDE

Col_VB=[1,length(Coord(1,:,1))];
Fil_VB=length(Coord(:,1,1));

ContFil=1;
 
for i=1:length(Col_VB)      %Marca Fila en matrix de coordenada Z
    for j=1:length(Coord(:,i))-1
        F_VBORDE(ContFil,1)=Coord(j,Col_VB(i),1);
        F_VBORDE(ContFil,2)=Coord(j,Col_VB(i),2);
        F_VBORDE(ContFil,3)=Coord(j,Col_VB(i),3);
        F_VBORDE(ContFil,4)=Coord(j+1,Col_VB(i),1);
        F_VBORDE(ContFil,5)=Coord(j+1,Col_VB(i),2);
        F_VBORDE(ContFil,6)=Coord(j+1,Col_VB(i),3);
        F_VBORDE(ContFil,7)=3;
        F_VBORDE(ContFil,8)=0;
        
        ContFil=ContFil+1;
    end
end

for i=1:length(Fil_VB)      %Marca Fila en matrix de coordenada Z
    for j=1:length(Coord(i,:,1))-1
        F_VBORDE(ContFil,1)=Coord(Fil_VB(i),j,1);
        F_VBORDE(ContFil,2)=Coord(Fil_VB(i),j,2);
        F_VBORDE(ContFil,3)=Coord(Fil_VB(i),j,3);
        F_VBORDE(ContFil,4)=Coord(Fil_VB(i),j+1,1);
        F_VBORDE(ContFil,5)=Coord(Fil_VB(i),j+1,2);
        F_VBORDE(ContFil,6)=Coord(Fil_VB(i),j+1,3);
        F_VBORDE(ContFil,7)=3;
        F_VBORDE(ContFil,8)=0;
        ContFil=ContFil+1;
    end
end


%% EXPORTADO DE DATOS

% CHEQUEOS
CTOTAL=M_PLACA+M_CARPETA+M_SISBER+M_EQUIPO+M_LUNIFC+M_VIGAS+M_CRAINE+M_MPILO;
MASA_TOTAL=[M_PLACA(:,1:3),(1/9.81).*CTOTAL(:,4:9)];
CHKMASA=(D_MPLACA+D_MCARPETA+D_MSISBER.*D_SISBER+D_MEQUIP+L_MUNIFC+D_MVIGAS+D_MCRAINE+MASA)./9.81;

CARGA_VERT=zeros(length(MASA_TOTAL(:,1)),length(MASA_TOTAL(1,:)));
CARGA_VERT(:,1:3)=MASA_TOTAL(:,1:3);
CARGA_VERT(:,6)=-9.81.*MASA_TOTAL(:,4);
CHKCV=(D_MPLACA+D_MCARPETA+D_MSISBER.*D_SISBER+D_MEQUIP+L_MUNIFC+D_MVIGAS+D_MCRAINE+MASA);

Ctot=sum(CTOTAL(:,4))

% ENCABEZADO GENERAL

% ENCABEZADO GENERAL

Enc0={  'PROYECTO:TESIS-ESCALAMIENTO_DE_SEÑALES_SISMICAS_EN_MUELLES';
        'AUTOR: ING. JUAN CARLOS PANTOJA';
        'ARCHIVO: MUELLE 20:1 - 323'};

%% EXPORTADO DE DATOS: ELEMENTOS
symbols =('A':'Z');

switch ctr2
   case 1
       
        % NODOS PARA APOYOS

        var=N_APOYO;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN NODOS: APOYOS PILOTES - LINKS'};
        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(N_APOYO(:,1))) ,' ELTOS');
            end
        end

        Enc(6,:)={'MÉTODO DE IMPORTACIÓN: '};
        Enc(6,:)=strcat(Enc(6,:),' PILOTES_','(7-',num2str(Exp_apoy(1)),')',',_LINK H1',' (',num2str(7+sum(Exp_apoy(1:2-1))),'-',num2str(Exp_apoy(2)),')',',_LINK H2',' (',num2str(7+sum(Exp_apoy(1:3-1))),'-',num2str(Exp_apoy(3)),')',', Cota Apoyos VC_ ',mat2str(Apoyos_Pilvc),', Cota Apoyos VG_ ',mat2str(Apoyos_Pilvg));

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('A. SUPPORTS - NODO.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS PARA LINKS H

        var=F_LINKH;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN NODOS: LINK H1 Y H2'};

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_LINKH(:,1))) ,' ELTOS');
            end
        end

        Enc(6,:)={'MÉTODO DE IMPORTACIÓN: '};
        for i=1:length(Exp_link)
            if i==1
                Enc(6,:)=strcat(Enc(6,:),Link_type(i),'_','(7-',num2str(Exp_link(i)),')');
            else
                Enc(6,:)=strcat(Enc(6,:),',_',Link_type(i),' (',num2str(7+sum(Exp_link(1:i-1))),'-',num2str(Exp_link(i)),')');
            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('B. LINK - SIMPLE BAR.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS PARA CONEXION

        var=F_CONEX;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN NODOS:'};
        for i=1:length(Link_type)+1
            if i==1
                Enc(4,:)=strcat(Enc(4,:),' Division Pile_ ',mat2str(PileDiv),', 3 CG x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');
            else 
                Enc(4,:)=strcat(Enc(4,:),num2str(Exp_link(i-1)),'_',Link_type(i-1),'x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');

            end
        end

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO ELEMENTOS CONEXION: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_CONEX(:,1))) ,' ELTOS');
            end
        end

        Enc(6,:)={'MÉTODO DE IMPORTACIÓN:'};
        for i=1:length(Exp_conex)
            if i==1
                Enc(6,:)=strcat(Enc(6,:),' PILOTE_ ',symbols(i),' (7-',num2str(Exp_conex(i)),')');
            else
                Enc(6,:)=strcat(Enc(6,:),', PILOTE_',symbols(i),' (',num2str(7+sum(Exp_conex(1:i-1))),'-',num2str(Exp_conex(i)),')');

            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('C. CONEXION - FRAME.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS PARA PILOTE

        var=F_PILE;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN NODOS:'};
        for i=1:length(Link_type)+1
            if i==1
                Enc(4,:)=strcat(Enc(4,:),' Division Pile_ ',mat2str(PileDiv),', 3 CG x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');
            else 
                Enc(4,:)=strcat(Enc(4,:),num2str(Exp_link(i-1)),'_',Link_type(i-1),'x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');

            end
        end

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO ELEMENTOS PILOTE: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_PILE(:,1))) ,' ELTOS');
            end
        end

        Enc(6,:)={'MÉTODO DE IMPORTACIÓN:'};
        for i=1:length(Exp_pile)
            if i==1
                Enc(6,:)=strcat(Enc(6,:),' PILOTE_ ',symbols(i),' (7-',num2str(Exp_pile(i)),')');
            else
                Enc(6,:)=strcat(Enc(6,:),', PILOTE_',symbols(i),' (',num2str(7+sum(Exp_pile(1:i-1))),'-',num2str(Exp_pile(i)),')');

            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('D1. PILE - FRAME.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS PARA PRESFUERZO

        var=F_PILE(:,1:6);    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN NODOS:'};
        for i=1:length(Link_type)+1
            if i==1
                Enc(4,:)=strcat(Enc(4,:),' 3 CG x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');
            else 
                Enc(4,:)=strcat(Enc(4,:),num2str(Exp_link(i-1)),'_',Link_type(i-1),'x_',num2str(length(MasaPilvc)+length(MasaPilvg)),' EJES;_');

            end
        end

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO ELEMENTOS PILOTE: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_PILE(:,1))) ,' ELTOS');
            end
        end

        Enc(6,:)={'MÉTODO DE IMPORTACIÓN:'};
        for i=1:length(Exp_pile)
            if i==1
                Enc(6,:)=strcat(Enc(6,:),' PILOTE_ ',symbols(i),' (7-',num2str(Exp_pile(i)),')');
            else
                Enc(6,:)=strcat(Enc(6,:),', PILOTE_',symbols(i),' (',num2str(7+sum(Exp_pile(1:i-1))),'-',num2str(Exp_pile(i)),')');

            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('D2. PRESSTRESS - FRAME .txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);


        % ELEMENTOS PARA FRAMES VIGA CABEZAL

        var=F_VCAB;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: FRAMES DE VIGAS CABEZAL'};

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_VCAB(:,1))) ,' ELTOS');
            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('E. VIGA CABEZAL - FRAME.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS PARA FRAMES VIGA GRUA

        var=F_VGRUA;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: FRAMES DE VIGAS GRUA'};

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(F_VGRUA(:,1))) ,' ELTOS');
            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('F. VIGA GRUA - FRAME.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);


        % ELEMENTOS PARA FRAMES VIGA DE BORDE

        var=F_VBORDE;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: FRAMES DE VIGAS DE BORDE'};

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                 Enc(5,:)=strcat( Enc(5,:),num2str(length(F_VBORDE(:,1))) ,' ELTOS');
            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('G. VIGA BORDE - FRAME.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % ELEMENTOS AREAS (SHELL)

        var=MatPuntos;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: ELEMENTOS TIPO SHELL - PLACA'};

        for i=1:2
            if i==1
                Enc(5,:)={'RECUENTO: '};
            else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(MatPuntos(:,1))) ,' ELTOS');
            end
        end

        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end
        fid = fopen('H. PLACA35CM - SHELL.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

    %% EXPORTADO DE CARGAS

        % MASA TOTAL
        
        var=MASA_TOTAL;    % VARIABLE A EXPORTAR
        

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: MASA TOTAL SIN MASA DE PILOTES'};
        
        for i=1:2
             if i==1
                Enc(5,:)={'RECUENTO: '};
             else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(MASA_TOTAL(:,1))) ,' ELTOS');
             end
        end
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('I. MASA_TOTAL - MASA.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);
        
        
        % CARGA VERTICAL D+0.1L 
        var=CARGA_VERT;    % VARIABLE A EXPORTAR

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: MASA TOTAL'};
        
        for i=1:2
             if i==1
                Enc(5,:)={'RECUENTO: '};
             else
                Enc(5,:)=strcat(Enc(5,:),num2str(length(CARGA_VERT(:,1))) ,' ELTOS');
             end
        end
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('J. CARGA D+0.1L - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        
    case 0
      
        % CARGAS: D_PLACA
        
        var=M_PLACA;    % VARIABLE A EXPORTAR
        
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_PLACA'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('K. D_PLACA - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        
        % CARGAS: D_CARPETA 
        
        var=M_CARPETA;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_CARPETA'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('L. D_CARPETA - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % CARGAS: D_SISTEMA BERTHING
        
        var=M_SISBER;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_SISTEMA BERTHING'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('L. D_SISTEMA BERTHING - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % CARGAS: D_EQUIPO
        
        var=M_EQUIPO;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        

        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_EQUIPO'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('M. D_EQUIP - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % CARGAS: L_UNIFORME_C
        
        var=M_LUNIFC;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS L_UNIFORME_C'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('N. L_UNIFORME_C - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);
        
        % CARGAS: D_VIGAS
        
        var=M_VIGAS;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_VIGAS'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('O. D_VIGAS - LOAD.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % CARGAS: D_CRAINE
        
        var=M_CRAINE;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: CARGAS D_CRAINE'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('P. D_CRAINE - NODAL MASS.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);

        % CARGAS: MASA PILOTES
        
        var=M_MPILO;    % VARIABLE A EXPORTAR
        var(:,6)=-1.*var(:,5);
        var(:,4:5)=0;
        
        Enc=Enc0;
        Enc(4,:)={'DESCRIPCIÓN: MASA DE NODOS EN PILOTES'};
        
        [nrows,ncols] = size(var);

        for i = 1:nrows
            for j=1:ncols
                if j==1
                    Enc(i+7,:)={num2str(var(i,1))};
                else
                    Enc(i+7,:)=strcat(Enc(i+7),',',num2str(var(i,j)));
                end
            end
        end

        fid = fopen('Q. MASA PILOTES - NODAL MASS.txt','w');
        fprintf(fid,'%s\r\n', Enc{:});
        fclose(fid);
       
end


        
        
            
        
        
