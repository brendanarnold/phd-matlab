function amp=nehint;

integ_tol=1e-6;

for z=1:10
    amp(z)=quad(@integrand,0,2*pi,integ_tol);
end

    function integrand_val=integrand(x);
        kf=1+0.1*sin(2*x).^2;
        vf=1+0.2*sin(2*x).^2;
        tau=(1+0.2*sin(2*x).^2)/(100+30*z.^2);
        integrand_val=kf.*vf.*tau;
    end
end