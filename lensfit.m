close all; clear all;
%% Load data
savefile = 'lensini.mat';
load(savefile);
l_OA = ER(1,1,3); % length of optical axis
%% extend lens to negative y plane
N_lens = extend_to_neg(N_lens);
X_lens = extend_to_neg(X_lens);


% update length of N
l_N = length(N_lens);
%% 2d to 3d
N_3d = zeros(3,l_N);
X_3d = zeros(3,l_N);
M_3d = zeros(3,l_N);
Y_3d = zeros(3,l_N);

N_3d(1:2,:,:) = N_lens;
X_3d(1:2,:,:) = X_lens;

M_3d(1,:,:) = N_lens(1,:,:);
M_3d(3,:,:) = N_lens(2,:,:);
Y_3d(1,:,:) = X_lens(1,:,:);
Y_3d(3,:,:) = X_lens(2,:,:);

%% Fit ellipse
mid_index = l_N/2;
EllipsParameters = zeros(mid_index,6);
for i = 1:mid_index
    n_s = N_3d(2:3,i);
    n_e = N_3d(2:3,mid_index+i);
    m_s = M_3d(2:3,i);
    m_e = M_3d(2:3,mid_index+i);
    extra = N_3d(2,i)*[cos(0.5),sin(0.5)]';
    XY = [m_s m_e n_s n_e extra]';
    EllipsParameters(i,:) = EllipseDirectFit(XY)';
end
%% plot
figure(1);
plot3(N_3d(1,:),N_3d(2,:),N_3d(3,:),'.r'); hold on;
plot3(X_3d(1,:),X_3d(2,:),X_3d(3,:),'.k');
plot3(M_3d(1,:),M_3d(2,:),M_3d(3,:),'.g');
plot3(Y_3d(1,:),Y_3d(2,:),Y_3d(3,:),'.k');
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-3 3]);
zlim([-3 3]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
grid on
hold off;

figure(2);
plot_points([m_s m_e n_s n_e],'.k');
xlim([-3 3]);
ylim([-3 3]);
xlabel('x [m]');
ylabel('y [m]');
grid on

figure(3);
hold on
for i = 1:length(EllipsParameters)
    %Convert the A to str 
    a = num2str(EllipsParameters(i,1)); 
    b = num2str(EllipsParameters(i,2));
    c = num2str(EllipsParameters(i,3));
    d = num2str(EllipsParameters(i,4));
    e = num2str(EllipsParameters(i,5));
    f = num2str(EllipsParameters(i,6));

    %Equation 
    figure(3);
    eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')']; 
    fimplicit(@(x,y) eval(eqt));
end
hold off
movegui
%% functions
function plot_point(P,c)
    plot(P(1),P(2),c)
end

function plot_points(P_vec,c)
    hold on;
    for P = P_vec       
        plot_point(P,c);
    end
    hold off;
end

function in = extend_to_neg(in)
    in2 = in;
    in2(2,:,:) = -in(2,:,:);
    in(:,:,(end+1):(2*length(in))) = in2;
end

function A = EllipseDirectFit(XY);
%
%  Direct ellipse fit, proposed in article
%    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%     "Direct Least Squares Fitting of Ellipses"
%     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%  Our code is based on a numerically stable version
%  of this fit published by R. Halir and J. Flusser
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%  This is a fast non-iterative ellipse fit.
%
%  It returns ellipses only, even if points are
%  better approximated by a hyperbola.
%  It is somewhat biased toward smaller ellipses.
%
centroid = mean(XY);   % the centroid of the data set
D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2];
D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,eval] = eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
A1 = evec(:,find(cond>0));
A = [A1; T*A1];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
     A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
A = A/norm(A);
end  %  EllipseDirectFit