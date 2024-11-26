% Define Proportional Controller
function u = PControl(r, y, us, Kc, umin, umax)
    e = r - y;
    v = us + Kc * e;
    u = max(umin, min(umax, v));
end