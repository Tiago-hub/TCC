function J = LQR_cal(x,U,Q,R,N,t)
    u=U.*ones(size(t));
    xt=x.'; %transposta de x
    ut=u.';
    for i=1:length(t)        
        a(i) = x(i,:)*Q*xt(:,i) + u(i,:)*R*ut(:,i) + 2*x(i,:)*N*ut(:,i);
        
    end
    J=trapz(a);

end