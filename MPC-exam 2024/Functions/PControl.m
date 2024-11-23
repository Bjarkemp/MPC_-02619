function [u] = PControl(r,z,us,Kc,umin,umax) 
        e = r-z; 
        v = us + Kc*e;
        u = max(umin,min(umax,v));
end