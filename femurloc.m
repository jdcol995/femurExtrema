function [ p1,p2 ] = femurloc( IMAGE )
%
 %
    %   femurloc es una funcion que toma como parametro una imgen 
    %   y le aplica un procesamiento para encontrar los extremos
    %   que delimitan la estrcutura osea del femur. La funcion regresa las
    %   coordenadas de los puntos de los extremos del femur. 
    %
    %
    %
    %   EJEMPLO:
    %      femur_img=dicomread('femur.dcm');
    %      [P1,P2]=femurloc(femur_img)
    %                            P1 =  55   176
    %                            P2 =  619  236     
    %
    %
    %
    %
warning off    
X_gray=rgb2gray(IMAGE);    
%%%%%%%%%%%%%%%%%%%%  Recorte Automatico   %%%%%%%%%%%%%%%%%%
    
DIM=size(X_gray);% se obtienen las dimensiones de la imagen
X=DIM(1)*0.1;    % 10% de las columnas derechas
Y=DIM(2)*0.15;   % 15% de los renglones superiores
X2=DIM(2)-DIM(2)*0.1;% 10% de las columnas izquierdas
   
%Se guarda la  ueva imagen en una matriz
img=X_gray(Y-1:DIM(1),X:X2);

%%%%%%%%%%%%%%%%    Fin del Recorte automatico    %%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
%%%%%%%%%%%%%%%           Pre procesamiento      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%aplicacion del filtro de difusion anisotropica
num_iter = 20;
delta_t = 1/4;
kappa = 30;
option = 2;
filt_img=anisodiff2D(img,num_iter,delta_t,kappa,option);
filt_img=uint8(filt_img);
%Se buscan los limites ideales para realizar el alargamiento al
%histograma
L_H=stretchlim(filt_img);
%Se ajusta la imagen a los limites encontrados en el histograma
im_adj=imadjust(filt_img,L_H);
 
%Se crea el elemnto estructurante tipo disco tama�o 25
SE=strel('disk',25);%Elemento estructurante
img_open=imopen(im_adj,SE);%Se realiza la apertura
  
%Se realiza la resta del la imagen original con la imagen abierta
imsubclos=imsubtract(im_adj,img_open);
   
% %Se crea un segundo elemento estrucurante tipo disco y tama�o 4
% SEE=strel('disk',4);%Segundo elemento estructurante
% %Se realiza el cierre de la imagen con el nuevo%SE
% imsubclos=imclose(imsub,SEE);
 
%Se obtiene el umbral optimo por metodo de Otsu
% Para este paso la imagen se divide en 4 secciones y se obtinen los
% umbrales de cada una, se elige el umbral de mayor valor como global
dim=size(imsubclos);
X=dim(1)*0.5;
Y=dim(2)*0.5;
crop1=imsubclos(1:X,1:Y);%       Primera seccion
crop2=imsubclos(1:X,Y-1:dim(2));%Segunda seccion    
crop3=imsubclos(X-1:dim(1),1:Y);%tercera seccion
crop4=imsubclos(X-1:dim(1),Y-1:dim(2));%cuarta seccion
 
T=graythresh(imsubclos);% Umbral de toda la imagen
%Umbrales de cada seccion
T1=graythresh(crop1);
T2=graythresh(crop2);
T3=graythresh(crop3);
T4=graythresh(crop4);
umbrales=[T,T1,T2,T3,T4];
%Eleccion del umbral de mayor tama�o
maximo=max(umbrales);
T=maximo;
%%%%%%%%%%%%%%%         Fin  Pre procesamiento      %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Se realiza la binarizacion de la imagen con el umbral de Otsu
bin=im2bw(imsubclos,T);
   
%Se realiza la apertura de imagen binaria utilizando un SE
%tipo disco tama�o 5
SE2=strel('square',9);
open_bin=imopen(bin,SE2);
%Se realiza el etiquetado de la imagen
labels=bwlabel(open_bin);
%        Se buscan las propiedades de las regiones etiquetadas
%intensidad media de las regiones
I_region = regionprops(labels, im_adj, 'MeanIntensity');
I_region=struct2cell(I_region);
I_region=cell2mat(I_region)';
%Centroide de las regiones etiquetadas
centroid_region = regionprops(labels,'centroid');%intensidad media de las regiones
centroid_region=cat(1, centroid_region.Centroid); 
 
%Se busca la cantidad maxima de etiquetas
labelMax=max(labels);
labelMax=max(labelMax);
 
%Vector para almacenar la distancia de todas las regiones
longitud=[];
%Vectores para almacenar las coordenadas de los puntos extremos
X1=[];%coordenadas x1
X2=[];%coordenadas x2
Y1=[];%coordenadas y1
Y2=[];%coordenadas y2
 
%        Vectores adicionales para la seleccion de la region
% Vector para almacenar las etiquetas de las regiones peque�as
small_Regions=[];
% Vector para almacenar las etiquetas de las regiones que son los
% suficientemente grande para ser femur
cand_regions=[];

%Longitud m�xima que se puede dar en la imagen
max_L=sqrt(((1-dim(2))^2)+((1-dim(1))^2));
  
%barrido a todas las etiquetas de la imagen
for i=1:labelMax
      %Se obtienen los renglones y columnas de la eqtiqueta 'i'
       [r,c]=find(labels==i);
       %coordenadas de cada punto de la region con etiqueta 'i'
       rc=[r c];
       %Se buscan las coordendas x1,y1 que son las coordenadas
       %de uno de los extremos del femur
       x1=rc(1,1);
       y1=rc(1,2);
       %las coordenadas del extremo de cada region son almacenadas
       %en dos vectores, uno para 'x' y otro para 'y'
       X1=[X1 x1];
       Y1=[Y1 y1];
       
       %Se buscan las coordendas del otro extremo de la region
       R=length(r);
       x2=rc(R,1);
       y2=rc(R,2);
       %Nuevamente, las coordendas del extrmo opuesto de la region
       %son almacencandas en dos vectores
       X2=[X2 x2];
       Y2=[Y2 y2];
       %Se mide la distancia entre lo extremos de la region y el resultado
       %se almacena en un vector
       dist=sqrt(((x1-x2)^2)+((y1-y2)^2));
       % Condicion para evaluar si la region es lo suficientemente
       % grande para ser considerada como candidato a femur
       if dist>(max_L*.21)
            longitud=[longitud dist];
            cand_regions=[cand_regions i];
       else
           small_Regions=[small_Regions i];
       end
end
   
% Se eliminan todos los indices correspondientes a las etiquetas de las
% regiones que son muy peque�as para ser candidatas a femur con lo que
% se quedan solo las etiquetas de las regiones grandes
X1(small_Regions)=[];
X2(small_Regions)=[];
Y1(small_Regions)=[];
Y2(small_Regions)=[];
I_region(small_Regions)=[];
centroid_region(small_Regions,:)=[];
for k=1: length(small_Regions)
        labels(labels==small_Regions(k))=0;
end

% Entropia de las regiones

entropia=[];
for k=1: length(cand_regions)
    aux_labels=labels;
    aux=cand_regions(k);
    aux_labels(aux_labels~=aux)=0;
    aux_labels=uint8(aux_labels);
    aux_labels=aux_labels/(max(max(aux_labels)));
    aux_im=aux_labels.*im_adj;

    p=imhist(aux_im(:));
    p(1)=[];
    p=p./sum(p);
    p(p==0)=[];
    H=-sum(p.*log2(p));
    entropia=[entropia H];
end
textura=1./entropia;

%        Seleccion de la region candidata a ser femur
scores=[];
for k=1: length(longitud)
       scores=[scores (1/4)*((I_region(k)/max(I_region))+(longitud(k)/max(longitud))+(centroid_region(k,2)/max(centroid_region(:,2))) + (textura(k)/max(textura)))];
end
femur_loc=find(scores==max(scores));

%         Union de regiones con caractaristicas similiares
% Limite de intesidad
intensity = I_region(femur_loc)-10;
% Limite de posicion 
centroid = abs(centroid_region(:,2)-centroid_region(femur_loc,2));
  
% Se buscan las regiones con caracteristicas
regions_i=find(I_region>intensity);
centroids_i=find(centroid<40);
   
% Si se encuentran caracterisitcas similares en ambos criterios la
% regiones se unen
if length(regions_i)==length(centroids_i)
  
     for k=1:length(regions_i)
         labels(labels==(cand_regions(regions_i(k))))=cand_regions(femur_loc);
     end
end

%        Adelgazamiento de la imagen
% Se eliminan todas las regiones que no pertenezcan a la que se
% considera femur
labels(labels~=cand_regions(femur_loc))=0;
% Se realiza el afelgazamiento
adelgado=bwmorph(labels,'thin',Inf);
% Se obtienen los puntos de la region adelgazada
[r,c]=find(adelgado==1);

% Se hace replica de las columnas y renglones de la region adelgazada
C=c;
R=r;
% Factor para submuestrar
DS=25;
% Submuestreo de la region adelgazada
C_D=downsample(C,DS);
R_D=downsample(R,DS);
% Interpolacion de la region submuestreada
xqn=min(C):0.1:max(C);
vq=interp1(C_D,R_D,xqn,'spline');
%     Division de la region en parte izquierda y derecha
% Parte izquierda
x_interI=xqn(1:round(length(xqn)/2));
y_interI=vq(1:round(length(vq)/2));

% Parte derecha
x_interD=xqn(round(length(xqn)/2):end);
y_interD=vq(round(length(vq)/2):end);

% Cambio de angulo cada 25 pixeles
intervalo=250;
anguloI=[];
anguloD=[];
new_XI=min(x_interI):max(x_interI);
new_XD=min(x_interD):max(x_interD);

for k=0:intervalo:length(x_interI)
    % Calculo del camio de angulo para la parte izquierda
    dy_I=y_interI(end)-y_interI(k+1);
    dx_I=x_interI(end)-x_interI(k+1);
    m_I=dy_I/dx_I;
    angulo=atan(m_I)*(180/pi);
    anguloI=[anguloI angulo];

    % Calculo del camio de angulo para la parte derecha
    if (intervalo+k)>length(x_interI)
        dy_D=y_interD(1)-y_interD(end);
        dx_D=x_interD(1)-x_interD(end);
        m_D=dy_D/dx_D;
        angulo=atan(m_D)*(180/pi);
        anguloD=[anguloD angulo];
        break
    end
    dy_D=y_interD(1)-y_interD(intervalo+k);
    dx_D=x_interD(1)-x_interD(intervalo+k);
    m_D=dy_D/dx_D;
    angulo=atan(m_D)*(180/pi);
    anguloD=[anguloD angulo];    
end

% Diferencia entre el cambio de angulo
diff_I=abs(diff(anguloI));
diff_D=abs(diff(anguloD));

% Umbral para realizar el corte
diff_I(diff_I<4)=0;
diff_D(diff_D<4)=0;

temp_I=diff_I(1:4);
temp_D=diff_D(end-3:end);

curv_I=find(temp_I);
curv_D=find(temp_D);

sz_I=length(curv_I);
sz_D=length(curv_D);

% Cantidad de pixeles a cortar en ambos extremos
pixels_RI=round(sz_I*(intervalo/10));
pixels_RD=round(sz_D*(intervalo/10));

% Coordenadas sobre las cuales marcar el punto
P1=[xqn(1)+pixels_RI round(vq(find(xqn==(xqn(1)+pixels_RI))))];
P2=[xqn(end)-pixels_RD round(vq(find(xqn==(xqn(end)-pixels_RD))))];

% Se grafica la imagen con los marcadores que denotan los extremos
% del femur
% Se guardan las poscisiones de las coordenadas de los extremos de la
% region adelgazada
% Compensasion de 10 pixeles
p1=[P1(1)+DIM(1)*0.1  P1(2)+DIM(2)*0.15];
p2=[P2(1)+DIM(1)*0.1  P2(2)+DIM(2)*0.15];
% p1=[P1(1)+DIM(1)*0.1  P1(2)+DIM(2)*0.15];
% p2=[P2(1)+DIM(1)*0.1  P2(2)+DIM(2)*0.15];

end

