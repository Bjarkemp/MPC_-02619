function [x_new] = qpsolver(H,g,low_x,up_x,mat,lb,ub,x_init)
    A_quad = [mat ; -mat];
    bound = [ub ; -lb];
    options=optimoptions("quadprog","Display","none");
    [x_new info] = quadprog(H,g,A_quad,bound,[],[],low_x,up_x,x_init,options);
end