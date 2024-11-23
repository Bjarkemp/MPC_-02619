function [u,i] = PIControl(i,r,z,us,Kc,taui,t_int,umin,umax) 
        e = r-z; 
        v = us + Kc*e + i;
        i = i+(Kc*t_int/taui)*e;
        u = max(umin,min(umax,v));
end