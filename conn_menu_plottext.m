function h=conn_menu_plottext(x,y,str,varargin)
% h=conn_menu_plottext(x,y,str)
% plots str along the coordinates of the curve (x,y)
% returns vector h of text handles (one for each letter/character in str)
%
%     % e.g.
%     fontsize=14;
%     str='En un lugar de la Mancha, de cuyo nombre no quiero acordarme, no ha mucho tiempo que vivía un hidalgo de los de lanza en astillero, adarga antigua, rocín flaco y galgo corredor. Una olla de algo más vaca que carnero, salpicón las más noches, duelos y quebrantos los sábados, lantejas los viernes, algún palomino de añadidura los domingos, consumían las tres partes de su hacienda. El resto della concluían sayo de velarte, calzas de velludo para las fiestas, con sus pantuflos de lo mesmo, y los días de entresemana se honraba con su vellorí de lo más fino. Tenía en su casa una ama que pasaba de los cuarenta y una sobrina que no llegaba a los veinte, y un mozo de campo y plaza que así ensillaba el rocín como tomaba la podadera. Frisaba la edad de nuestro hidalgo con los cincuenta años. Era de complexión recia, seco de carnes, enjuto de rostro, gran madrugador y amigo de la caza. Quieren decir que tenía el sobrenombre de «Quijada», o «Quesada», que en esto hay alguna diferencia en los autores que deste caso escriben, aunque por conjeturas verisímiles se deja entender que se llamaba «Quijana». Pero esto importa poco a nuestro cuento: basta que en la narración dél no se salga un punto de la verdad.';
%     clf;
%     axis equal tight;
%     set(gca,'xlim',[-6 17],'ylim',[-6 6]);
%     p=5*linspace(0,1,1e3);
%     x=cos(-2*pi*p).*(.5+p);
%     y=sin(-2*pi*p).*(.5+p);
%     x=[x 11-flip(x)];
%     y=[y -flip(y)];
%     h=conn_menu_plottext(x,y,str,'fontsize',fontsize);    

va='bottom'; 
if nargin<1||isempty(x), x=0; end
if nargin<2||isempty(y), y=0; end
if nargin<3||isempty(str),str=''; end
if numel(x)==1&&numel(y)==1,
    a=-2*pi*(-49.5:49.5)/100+angle(x+1i*y); % note: single coordinates interpreted as circle around zero
    b=abs(x+1i*y);
    x=b*cos(a);
    y=b*sin(a);
    if x(1)>x(end), x=fliplr(x); y=fliplr(y); va='top'; end
elseif numel(x)==1&&numel(y)>1, x=x+zeros(size(y));  
elseif numel(x)>1&&numel(y)==1, y=y+zeros(size(x));  
end

letterextent=[0.2069 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2069 0.2414 0.2414 0.2759 0.2759 0.2414 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2759 0.2069 0.2759 0.2759 0.2414 0.2414 0.3793 0.4828 0.4828 0.8276 0.5517 0.2414 0.2414 0.2414 0.3103 0.5172 0.2414 0.3448 0.2414 0.2759 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.2414 0.2414 0.5172 0.5172 0.5172 0.4828 0.6897 0.5517 0.5862 0.6207 0.5862 0.5172 0.4828 0.6552 0.6207 0.2414 0.4483 0.5862 0.4828 0.7241 0.6207 0.6552 0.5517 0.6552 0.5862 0.5517 0.4828 0.6207 0.5172 0.7931 0.5172 0.5517 0.5172 0.2414 0.2759 0.2414 0.5172 0.4138 0.2069 0.4483 0.5172 0.4483 0.5172 0.4483 0.2759 0.4828 0.4828 0.2069 0.2069 0.4483 0.2069 0.7241 0.4828 0.4828 0.5172 0.5172 0.2759 0.4138 0.2759 0.4828 0.4138 0.6552 0.4483 0.4138 0.4138 0.2759 0.2069 0.2759 0.5172 0.3448 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.4138 0.2414 0.2414 0.4828 0.4828 0.4828 0.4828 0.2069 0.4828 0.2069 0.6897 0.3448 0.4138 0.5172 0.3448 0.6897 0.2069 0.3448 0.5172 0.2759 0.2759 0.2069 0.4828 0.5172 0.2414 0.2069 0.2759 0.3448 0.4138 0.7241 0.7241 0.7241 0.4828 0.5517 0.5517 0.5517 0.5517 0.5517 0.5517 0.7931 0.6207 0.5172 0.5172 0.5172 0.5172 0.2414 0.2414 0.2414 0.2414 0.5862 0.6207 0.6552 0.6552 0.6552 0.6552 0.6552 0.5172 0.6552 0.6207 0.6207 0.6207 0.6207 0.5517 0.5517 0.4828 0.4483 0.4483 0.4483 0.4483 0.4483 0.4483 0.7241 0.4483 0.4483 0.4483 0.4483 0.4483 0.2069 0.2069 0.2069 0.2069 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.4828 0.5172 0.4828 0.4828 0.4828 0.4828 0.4828 0.4138 0.5172 0.4138];
haxis=find(cellfun(@ischar,varargin)); haxis=haxis(find(strcmpi(varargin(haxis),'parent'),1));
if isempty(haxis), haxis=gca; 
else haxis=varargin{haxis+1};
end
hold(haxis,'on');
halign=find(cellfun(@ischar,varargin)); halign=halign(find(strcmpi(varargin(halign),'horizontalalignment'),1));
if isempty(halign), halign='left';
else halign=lower(varargin{halign+1});
end
h0=text(x(1),y(1),str,varargin{:},'parent',haxis);
h1=plot(x,y,'.','parent',haxis);
%hold(haxis,'off');
ext=get(h0,'extent');
L0=ext(3);
L1=ext(4);
dL=sqrt(max(0,reshape(x(2:end)-x(1:end-1),1,[]).^2+reshape(y(2:end)-y(1:end-1),1,[]).^2));
dL=[0 dL];
L=cumsum(dL);
maxL=L(end);
letterwidth=letterextent(max(1,min(numel(letterextent),double(str))));
slw=sum(letterwidth);
mlw=mean(letterwidth);
%L0=min(maxL,L0);
letterposition=L0*(cumsum(letterwidth)-letterwidth/2)/slw;
switch(halign),
    case 'left',
    case 'center', letterposition=letterposition-mean(letterposition)+interp1(1:numel(L),L,(1+numel(L))/2,'linear');
    case 'right',  letterposition=letterposition-max(letterposition)+L(end);
end
position_x = interp1(L,x,letterposition,'linear','extrap');
position_y = interp1(L,y,letterposition,'linear','extrap');
dl=L0*letterwidth/2/slw;
rotation = 180/pi*(0+angle(interp1(L,x,letterposition+dl,'linear','extrap')-interp1(L,x,letterposition-dl,'linear','extrap') + 1i*(interp1(L,y,letterposition+dl,'linear','extrap')-interp1(L,y,letterposition-dl,'linear','extrap'))));
h=[];
for n=1:numel(str)
    h(n)=text(position_x(n),position_y(n),str(n),varargin{:},'horizontalalignment','center','verticalalignment',va,'interpreter','none','rotation',rotation(n),'parent',haxis);
end
delete([h0 h1]);
