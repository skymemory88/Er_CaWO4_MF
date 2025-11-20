function [Jx, Jy, Jz, Ix, Iy, Iz, Jxh, Jyh, Jzh, Ixh, Iyh, Izh] = spin_operators(J, I)
    % Electronic spin operators
    Jz = diag(J:-1:-J);
    Jp = diag(sqrt((J - (J-1:-1:-J)).*(J+1 + (J-1:-1:-J))), 1);
    Jm = Jp';
    Jx = (Jp+Jm)/2;
    Jy = (Jp-Jm)/2i;
    
    % Nuclear spin operators
    Iz = diag(I:-1:-I);
    Ip = diag(sqrt((I - (I-1:-1:-I)) .* (I+1 + (I-1:-1:-I))), 1);
    Im = Ip';
    Ix = (Ip + Im)/2;
    Iy = (Ip - Im)/2i;
    
    % Hybrid operators
    Jxh = kron(Jx,eye(2*I+1));
    Jyh = kron(Jy,eye(2*I+1));
    Jzh = kron(Jz,eye(2*I+1));
    
    Ixh = kron(eye(2*J+1),Ix);
    Iyh = kron(eye(2*J+1),Iy);
    Izh = kron(eye(2*J+1),Iz);
end