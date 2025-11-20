function [Hcf] = cf(J, B, rot)
rot = rot*pi/180; % degree to rad
rot_opt = 'coefficient'; % rotation by 'operator'/'coefficient'
%----- Form the J operators
n = 2*J+1;
Jz = diag(J:-1:-J);
Jp = diag(sqrt((J-(J-1:-1:-J)).*(J+1+(J-1:-1:-J))),1);
Jm = Jp';
% Jx = (Jp+Jm)/2;
% Jy = (Jp-Jm)/2i;

%----- Build Vcf from the stevens operators:
X = J*(J+1);
O20 = 3*Jz^2 - X*eye(n);

O40 = 35*Jz^4-(30*X-25)*Jz^2 + (3*X^2-6*X)*eye(n);
O44c = (Jp^4+Jm^4)/2;
O44s = (Jp^4-Jm^4)/(1i*2);

O60 = 231*Jz^6 - (315*X-735)*Jz^4 + (105*X^2-525*X+294)*Jz^2 + (-5*X^3+40*X^2-60*X)*eye(n);
O64c = ((11*Jz^2-X*eye(n)-38*eye(n)) * (Jp^4+Jm^4) + (Jp^4+Jm^4) * (11*Jz^2-X*eye(n)-38*eye(n)))/4;
O64s = ((11*Jz^2-X*eye(n)-38*eye(n)) * (Jp^4-Jm^4) + (Jp^4-Jm^4) * (11*Jz^2-X*eye(n)-38*eye(n)))/(1i*4);  

% Crystal field rotation
switch rot_opt
    case 'operator' % rotate crystal field (method 1)
        Ur = [cos(4*rot)  -sin(4*rot)
              sin(4*rot)  +cos(4*rot)];
        Ur = kron(Ur,eye(2*J+1));
        O44 = [O44c O44s]';
        O44 = (Ur*O44)'; % method 1: rotate Operators
        O44c = O44(:,1:2*J+1);
        O44s = O44(:,2*J+2:end);
        O64 = [O64c O64s]';
        O64 = (Ur*O64)'; % method 1: rotate Operators
        O64c = O64(:,1:2*J+1);
        O64s = O64(:,2*J+2:end);
        
        % alternative: explicit form
        % O44co = O44c; O44so = O44s;
        % O64co = O64c; O64so = O64s;
        % O44c = cos(4*rot)*O44co - sin(4*rot)*O44so;
        % O44s = sin(4*rot)*O44co + cos(4*rot)*O44so;
        % O64c = cos(4*rot)*O64co - sin(4*rot)*O64so;
        % O64s = sin(4*rot)*O64co + cos(4*rot)*O64so;
    case 'coefficient' % rotate crystal field (method 2)
        Br = [ cos(4*rot)  sin(4*rot)
              -sin(4*rot)  cos(4*rot)];
        B44 = Br * ([B(3) B(4)]');
        B44 = B44';
        B(3) = B44(1);
        B(4) = B44(2);
        
        B64 = Br * ([B(6) B(7)]');
        B64 = B64';
        B(6) = B64(1);
        B(7) = B64(2);
        
        % alternative: explicit form
        % B44c = B(3); B44s = B(4);
        % B64c = B(6); B64s = B(7);
        % B(3) = cos(4*rot)*B44c + sin(4*rot)*B44s;
        % B(4) = -sin(4*rot)*B44c + cos(4*rot)*B44s;
        % B(6) = cos(4*rot)*B64c + sin(4*rot)*B64s;
        % B(7) = -sin(4*rot)*B64c + cos(4*rot)*B64s;
    otherwise
        % No rotation
end
Hcf = B(1)*O20 + B(2)*O40 + B(3)*O44c + B(4)*O44s + B(5)*O60 + B(6)*O64c + B(7)*O64s;
end