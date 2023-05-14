function [tt, r, s, tt_adf, r_adf, s_adf] = DR(x,y,z)
[T,N,B,k,t] = frenet(x,y,z);
% plot3(x,y,z), hold on
% quiver3(x',y',z',T(:,1),T(:,2),T(:,3),'color','r')
% quiver3(x',y',z',N(:,1),N(:,2),N(:,3),'color','g')
% quiver3(x',y',z',B(:,1),B(:,2),B(:,3),'color','b')
% legend('Curve','Tangent','Normal','Binormal')

n = size(x,2);
boundarydata0 = [x(:,1) y(:,1) z(:,1)];
boundarydata1 = [x(:,n) y(:,n) z(:,n)];
% U0 = [B(1,:);N(1,:);T(1,:)];

% doubleReflection
[t1, r1, s1] = Frame(boundarydata0, T(1,:));
XX = [x' y' z'];
tt1 = t1;
rr1 = r1;
ss1 = s1;
for i=1:n-1
    v1 = XX(i+1,:) - XX(i,:);
    c1 = dot(v1,v1);
    ref_L_i = rr1(i,:) - (2/c1) * dot(v1,rr1(i,:)) * v1;
    tan_L_i = tt1(i,:) - (2/c1) * dot(v1,tt1(i,:)) * v1;
    v2 = T(i+1,:) - tan_L_i;
    c2 = dot(v2,v2);
    ref_next = ref_L_i - (2 / c2) * dot(v2, ref_L_i) * v2;
    [tt_n, r_n, s_n] = Frame(ref_next,T(i+1,:));
    tt = [tt1;tt_n];
    tt1 = tt;
    r = [rr1;r_n];
    rr1 = r;
    s = [ss1;s_n];
    ss1 = s;
end
% figure
% plot3(x,y,z), hold on
% quiver3(x',y',z',tt(:,1),tt(:,2),tt(:,3),'color','r')
% quiver3(x',y',z',r(:,1),r(:,2),r(:,3),'color','g')
% quiver3(x',y',z',s(:,1),s(:,2),s(:,3),'color','b')
% legend('Curve','Tangent','Normal','Binormal')

% adjustFrames
[t_last, r_last, s_last] = Frame(boundarydata1, T(n,:));
maxAngle = -diff_angle(r(n,:), r_last);
parameterValues = 0;
step = 1 / (n-1);
tt_adf1 = [];
r_adf1 = [];
s_adf1 = [];
for j=1:n
    adf = adjustmentFunc(parameterValues, maxAngle, 0);  %twistCount=0
    parameterValues = parameterValues + step;
    ro = rotate_my(r(j,:), adf, tt(j,:));
    [tt_adf, r_adf, s_adf] = Frame(ro,tt(j,:));
    tt_adf = [tt_adf1;tt_adf];
    tt_adf1 = tt_adf;
    r_adf = [r_adf1;r_adf];
    r_adf1 = r_adf;
    s_adf = [s_adf1;s_adf];
    s_adf1 = s_adf;
end
% figure
% plot3(x,y,z), hold on
% quiver3(x',y',z',tt_adf(:,1),tt_adf(:,2),tt_adf(:,3),'color','r')
% quiver3(x',y',z',r_adf(:,1),r_adf(:,2),r_adf(:,3),'color','g')
% quiver3(x',y',z',s_adf(:,1),s_adf(:,2),s_adf(:,3),'color','b')
% legend('Curve','Tangent','Normal','Binormal')
end


function [tt, r, s] = Frame(reference,tangent)
tt = norm_my(tangent);
p_r_t = dot(reference,tangent)./dot(tangent,tangent).*tangent;
r = norm_my(reference - p_r_t);
s = norm_my(cross(tt, r));
end

function T = norm_my(V)
temp = sqrt( sum(V.^2,2) );
T = V ./ repmat(temp, 1, 3);
end

function b = diff_angle(a, other)
aa = norm_my(a);
oo = norm_my(other);
a = dot(aa,oo);
if a>1
    b = 0;
elseif a < -1
    b = pi;
else
    b = acos(a);
end     
end

function o = adjustmentFunc(t, maxAngle, twistCount)
o = t *(maxAngle + twistCount*2*pi);
end

function vector = rotate_my(r, angle, axis)
ea = exist('axis');
if ea == 0
    u = [0 0 1];
else
    u = norm_my(axis);
end
c = cos(angle);
s = sin(angle);
t = 1-c;
x = u(1);
y = u(2);
z = u(3);
m11=t*x*x+c;
m12=t*x*y-z*s;
m13=t*x*z+y*s;
m21=t*x*y+z*s;
m22=t*y*y+c;
m23=t*y*z-x*s;
m31=t*x*z-y*s;
m32=t*y*z+x*s;
m33=t*z*z+c;
sx=r(1);
sy=r(2);
sz=r(3);
vector = [(m11*sx+m12*sy+m13*sz), (m21*sx+m22*sy+m23*sz), (m31*sx+m32*sy+m33*sz)];
end