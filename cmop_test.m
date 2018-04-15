
function fobj = cmop_test(name)

switch name
    case 'LIRCMOP1'
        fobj = @LIRCMOP1;
    case 'LIRCMOP2'
        fobj = @LIRCMOP2;
    case 'LIRCMOP3'
        fobj = @LIRCMOP3;
    case 'LIRCMOP4'
        fobj = @LIRCMOP4;
    case 'LIRCMOP5'
        fobj = @LIRCMOP5;
    case 'LIRCMOP6'
        fobj = @LIRCMOP6;
    case 'LIRCMOP7'
        fobj = @LIRCMOP7;
    case 'LIRCMOP8'
        fobj = @LIRCMOP8;
    case 'LIRCMOP9'
        fobj = @LIRCMOP9;
    case 'LIRCMOP10'
        fobj = @LIRCMOP10;
    case 'LIRCMOP11'
        fobj = @LIRCMOP11;
    case 'LIRCMOP12'
        fobj = @LIRCMOP12;
    case 'LIRCMOP13'
        fobj = @LIRCMOP13;
    case 'LIRCMOP14'
        fobj = @LIRCMOP14;
    case 'CF1'
        fobj = @CF1;
    case 'CF2'
        fobj = @CF2;
    case 'CF3'
        fobj = @CF3;
    case 'CF4'
        fobj = @CF4;
    case 'CF5'
        fobj = @CF5;
    otherwise
        error('The optimized problem is not exsit!');
end
end


%% LIRCMOP1
function [y,c] = LIRCMOP1(x)
x_odd = x(3:2:end,:); x_even = x(2:2:end,:);
len_odd = size(x_odd,1); len_even = size(x_even,1);
g_1 = sum((x_odd - repmat(x(1,:),len_odd,1)).^2,1);
g_2 = sum((x_even - repmat(x(1,:),len_even,1)).^2,1);

y(1,:) = x(1,:) + g_1;
y(2,:) = 1 - x(1,:) .^ 2 + g_2;

% Constrains
c(1,:) = (g_1 - 0.5).*(0.51 - g_1);
c(2,:) = (g_2 - 0.5).*(0.51 - g_2);

end

%% LIRCMOP2
function [y,c] = LIRCMOP2(x)
x_odd = x(3:2:end,:); x_even = x(2:2:end,:);
len_odd = size(x_odd,1); len_even = size(x_even,1);
g_1 = sum((x_odd - repmat(x(1,:),len_odd,1)).^2,1);
g_2 = sum((x_even - repmat(x(1,:),len_even,1)).^2,1);

y(1,:) = x(1,:) + g_1;
y(2,:) = 1 - sqrt(x(1,:)) + g_2;

% Constrains
c(1,:) = (g_1 - 0.5).*(0.51 - g_1);
c(2,:) = (g_2 - 0.5).*(0.51 - g_2);

end

%% LIRCMOP3
function [y,c] = LIRCMOP3(x)
x_odd = x(3:2:end,:); x_even = x(2:2:end,:);
len_odd = size(x_odd,1); len_even = size(x_even,1);
g_1 = sum((x_odd - repmat(x(1,:),len_odd,1)).^2,1);
g_2 = sum((x_even - repmat(x(1,:),len_even,1)).^2,1);

y(1,:) = x(1,:) + g_1;
y(2,:) = 1 - x(1,:) .^ 2 + g_2;

% Constrains
c(1,:) = (g_1 - 0.5).*(0.51 - g_1);
c(2,:) = (g_2 - 0.5).*(0.51 - g_2);
c(3,:) = sin(20 * pi * x(1,:)) - 0.5;

end

%% LIRCMOP4
function [y,c] = LIRCMOP4(x)
x_odd = x(3:2:end,:); x_even = x(2:2:end,:);
len_odd = size(x_odd,1); len_even = size(x_even,1);
g_1 = sum((x_odd - repmat(x(1,:),len_odd,1)).^2,1);
g_2 = sum((x_even - repmat(x(1,:),len_even,1)).^2,1);

y(1,:) = x(1,:) + g_1;
y(2,:) = 1 - sqrt(x(1,:)) + g_2;

% Constrains
c(1,:) = (g_1 - 0.5).*(0.51 - g_1);
c(2,:) = (g_2 - 0.5).*(0.51 - g_2);
c(3,:) = sin(20 * pi * x(1,:)) - 0.5;

end


%% LIRCMOP5
function [y,c] = LIRCMOP5(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end
gx=0.7057;

y(1,:)=x(1,:)+10*sum1+gx;
y(2,:)=1-x(1,:).^0.5+10.*sum2+gx;

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=[1.6,2.5];
q=[1.6,2.5];
a=[2,2];
b=[4,8];
r=0.1;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
for k=1:2
    c(k,:)=((y(1,:)-p(k))*cos(theta)-(y(2,:)-q(k))*sin(theta)).^2/(a(k)^2)+...
        ((y(1,:)-p(k))*sin(theta)+(y(2,:)-q(k))*cos(theta)).^2/(b(k)^2)-r;
    
end
end

%% LIRCMOP6
function [y,c] = LIRCMOP6(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end
gx=0.7057;

y(1,:)=x(1,:)+10*sum1+gx;
y(2,:)=1-x(1,:).^2+10.*sum2+gx;

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=[1.8,2.8];
q=[1.8,2.8];
a=[2,2];
b=[8,8];
r=0.1;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
for k=1:2
    c(k,:)=((y(1,:)-p(k))*cos(theta)-(y(2,:)-q(k))*sin(theta)).^2/(a(k)^2)+...
        ((y(1,:)-p(k))*sin(theta)+(y(2,:)-q(k))*cos(theta)).^2/(b(k)^2)-r;
    
end
end

%% LIRCMOP7
function [y,c] = LIRCMOP7(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end
gx=0.7057;

y(1,:)=x(1,:)+10*sum1+gx;
y(2,:)=1-x(1,:).^0.5+10.*sum2+gx;

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=[1.2,2.25,3.5];
q=[1.2,2.25,3.5];
a=[2,2.5,2.5];
b=[6,12,10];
r=0.1;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
for k=1:3
    c(k,:)=((y(1,:)-p(k))*cos(theta)-(y(2,:)-q(k))*sin(theta)).^2/(a(k)^2)+...
        ((y(1,:)-p(k))*sin(theta)+(y(2,:)-q(k))*cos(theta)).^2/(b(k)^2)-r;
    
end
end

%% LIRCMOP8
function [y,c] = LIRCMOP8(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end
gx=0.7057;

y(1,:)=x(1,:)+10*sum1+gx;
y(2,:)=1-x(1,:).^2+10.*sum2+gx;

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=[1.2,2.25,3.5];
q=[1.2,2.25,3.5];
a=[2,2.5,2.5];
b=[6,12,10];
r=0.1;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
for k=1:3
    c(k,:)=((y(1,:)-p(k))*cos(theta)-(y(2,:)-q(k))*sin(theta)).^2/(a(k)^2)+...
        ((y(1,:)-p(k))*sin(theta)+(y(2,:)-q(k))*cos(theta)).^2/(b(k)^2)-r;
    
end
end

%% LIRCMOP9
function [y,c] = LIRCMOP9(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end

y(1,:)=1.7057*x(1,:).*(10*sum1+1);
y(2,:)=1.7057*(1-x(1,:).^2).*(10*sum2+1);

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=1.4;
q=1.4;
a=1.5;
b=6;
r=0.1;alpha=0.25*pi;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
c(1,:)=(((y(1,:)-p)*cos(theta)-(y(2,:)-q)*sin(theta)).^2)/(a^2)+...
    (((y(1,:)-p)*sin(theta)+(y(2,:)-q)*cos(theta)).^2)/(b^2)-r;
c(2,:)=y(1,:)*sin(alpha)+y(2,:)*cos(alpha)-sin(4*pi*(y(1,:)*cos(alpha)-y(2,:)*sin(alpha)))-2;


end

%% LIRCMOP10
function [y,c] = LIRCMOP10(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end

y(1,:)=1.7057*x(1,:).*(10*sum1+1);
y(2,:)=1.7057*(1-x(1,:).^0.5).*(10*sum2+1);

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=1.1;
q=1.2;
a=2;
b=4;
r=0.1;alpha=0.25*pi;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
c(1,:)=(((y(1,:)-p)*cos(theta)-(y(2,:)-q)*sin(theta)).^2)/(a^2)+...
    (((y(1,:)-p)*sin(theta)+(y(2,:)-q)*cos(theta)).^2)/(b^2)-r;
c(2,:)=y(1,:)*sin(alpha)+y(2,:)*cos(alpha)-sin(4*pi*(y(1,:)*cos(alpha)-y(2,:)*sin(alpha)))-1;


end

%% LIRCMOP11
function [y,c] = LIRCMOP11(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end

y(1,:)=1.7057*x(1,:).*(10*sum1+1);
y(2,:)=1.7057*(1-x(1,:).^0.5).*(10*sum2+1);

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=1.2;
q=1.2;
a=1.5;
b=5;
r=0.1;alpha=0.25*pi;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
c(1,:)=(((y(1,:)-p)*cos(theta)-(y(2,:)-q)*sin(theta)).^2)/(a^2)+...
    (((y(1,:)-p)*sin(theta)+(y(2,:)-q)*cos(theta)).^2)/(b^2)-r;
c(2,:)=y(1,:)*sin(alpha)+y(2,:)*cos(alpha)-sin(4*pi*(y(1,:)*cos(alpha)-y(2,:)*sin(alpha)))-2.1;


end

%% LIRCMOP12
function [y,c] = LIRCMOP12(x)
variable_length=size(x,1);
popsize=size(x,2);
sum1=zeros(1,popsize);
sum2=zeros(1,popsize);
for j=2:variable_length
    if mod(j,2)==1
        sum1=sum1+(x(j,:)-sin((0.5*j/variable_length*pi)*x(1,:))).^2;
    else
        sum2=sum2+(x(j,:)-cos((0.5*j/variable_length*pi)*x(1,:))).^2;
    end
end

y(1,:)=1.7057*x(1,:).*(10*sum1+1);
y(2,:)=1.7057*(1-x(1,:).^2).*(10*sum2+1);

%%%%%%%%% constraints %%%%%%%%%%%%%%%
%%%%%%%%% parameters %%%%%%%%%%%%%%%%
p=1.6;
q=1.6;
a=1.5;
b=6;
r=0.1;alpha=0.25*pi;theta=-0.25*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(2,popsize);
c(1,:)=(((y(1,:)-p)*cos(theta)-(y(2,:)-q)*sin(theta)).^2)/(a^2)+...
    (((y(1,:)-p)*sin(theta)+(y(2,:)-q)*cos(theta)).^2)/(b^2)-r;
c(2,:)=y(1,:)*sin(alpha)+y(2,:)*cos(alpha)-sin(4*pi*(y(1,:)*cos(alpha)-y(2,:)*sin(alpha)))-2.5;


end

%% LIRCMOP13
function [y,c] = LIRCMOP13(x)
variable_length=size(x,1);
popsize=size(x,2);
sum=zeros(1,popsize);
for j=3:variable_length
    sum=sum+10*(x(j,:)-0.5).^2;
end

y(1,:)=(1.7057+sum).*cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:));
y(2,:)=(1.7057+sum).*cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:));
y(3,:)=(1.7057+sum).*sin(0.5*pi*x(1,:));
%%%%%%%%% constraints %%%%%%%%%%%%%%%
gx=y(1,:).^2+y(2,:).^2+y(3,:).^2;
c=zeros(2,popsize);
c(1,:)=(gx-9).*(gx-4);
c(2,:)=(gx-3.61).*(gx-3.24);

end

%% LIRCMOP14
function [y,c] = LIRCMOP14(x)
variable_length=size(x,1);
popsize=size(x,2);
sum=zeros(1,popsize);
for j=3:variable_length
    sum=sum+10*(x(j,:)-0.5).^2;
end

y(1,:)=(1.7057+sum).*cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:));
y(2,:)=(1.7057+sum).*cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:));
y(3,:)=(1.7057+sum).*sin(0.5*pi*x(1,:));
%%%%%%%%% constraints %%%%%%%%%%%%%%%
gx=y(1,:).^2+y(2,:).^2+y(3,:).^2;
c=zeros(2,popsize);
c(1,:)=(gx-9).*(gx-4);
c(2,:)=(gx-3.61).*(gx-3.24);
c(3,:)=(gx-3.0625).*(gx-2.56);
end

%% CF1
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF1(x)
    a            = 1.0;
    N            = 10.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:) + 2.0*tmp2/size(2:2:dim,2);
    c(1,:)       = y(1,:) + y(2,:) - a*abs(sin(N*pi*(y(1,:)-y(2,:)+1.0))) - 1.0;
    clear Y;
end

%% CF2
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF2(x)
    a           = 1.0;
    N           = 2.0;
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
    t           = y(2,:) + sqrt(y(1,:)) - a*sin(N*pi*(sqrt(y(1,:))-y(2,:)+1.0)) - 1.0;
    c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
    clear tmp;
end

%% CF3
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF3(x)
    a            = 1.0;
    N            = 2.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1         = zeros(dim,num);
    tmp1(2:dim,:)= Y(2:dim,:).^2;
    tmp2         = zeros(dim,num);
    tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
    tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
    tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
    y(1,:)       = x(1,:)          + 2.0*tmp11/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:).^2 + 2.0*tmp21/size(2:2:dim,2);
    c(1,:)       = y(2,:) + y(1,:).^2 - a*sin(N*pi*(y(1,:).^2-y(2,:)+1.0)) - 1.0;   
    clear Y tmp1 tmp2;
end

%% CF4
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF4(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
    tmp2        = sum(tmp(4:2:dim,:).^2);  % even index
    index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
    index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
    tmp(2,index1) = abs(tmp(2,index1));
    tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
    y(1,:)      = x(1,:)                  + tmp1;
    y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
    t           = x(2,:) - sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
    c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
    clear tmp index1 index2;
end

%% CF5
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF5(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
    tmp2        = sum(2.0*tmp(4:2:dim,:).^2-cos(4.0*pi*tmp(4:2:dim,:))+1.0);  % even index
    index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
    index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
    tmp(2,index1) = abs(tmp(2,index1));
    tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
    y(1,:)      = x(1,:)                  + tmp1;
    y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
    c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
    clear tmp;
end

%% CF6
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF6(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
    tmp2        = sum(tmp(2:2:dim,:).^2);  % even index
    y(1,:)      = x(1,:)            + tmp1;
    y(2,:)      = (1.0 - x(1,:)).^2 + tmp2;
    tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
    c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
    tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
    c(2,:)      = x(4,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
    clear tmp;
end

%% CF7
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF7(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp2        = sum(2.0*tmp(6:2:dim,:).^2-cos(4.0*pi*tmp(6:2:dim,:))+1.0);  % even index
    tmp(2,:)    = tmp(2,:).^2;
    tmp(4,:)    = tmp(4,:).^2;
    y(1,:)      = x(1,:)                                  + tmp1;
    y(2,:)      = (1.0 - x(1,:)).^2 + tmp(2,:) + tmp(4,:) + tmp2;
    tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
    c(1,:)      = x(2,:) - sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
    tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
    c(2,:)      = x(4,:) - sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
    clear tmp;
end

%% CF8
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF8(x)
    N           = 2.0;
    a           = 4.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*abs(sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0))) - 1.0;
    clear Y;
end

%% CF9
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF9(x)
    N           = 2.0;
    a           = 3.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
    clear Y;
end

%% CF10
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF10(x)
    a           = 1.0;
    N           = 2.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]));
    H           = zeros(dim,num);
    H(3:dim,:)  = 4.0*Y(3:dim,:).^2 - cos(8.0*pi*Y(3:dim,:)) + 1.0;
    tmp1        = sum(H(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(H(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(H(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
    clear Y H;
end
