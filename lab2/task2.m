close all
mu = 6;
beta = 0.01;
x0 = [900;100];
tot_pop = sum(x0);
t_start = 0;
t_end = 20;
dt = 0.0001;

ode = @(t,x) [
    mu*x(2) - beta*x(2)*x(1);
    -mu*x(2) + beta*x(1)*x(2);
];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4, 'Stats', 'on');

[t_stoch,x_stoch] = stoch_sim(x0,mu,beta,t_start,t_end,dt);
[t,x] = ode45(ode, [t_start t_end], x0, options);
plot(t_stoch,x_stoch./sum(x0));
hold on
plot(t,x./sum(x0));
legend("\rho_s stochastic", "\rho_i stochastic", "\rho_s", "\rho_i", "FontSize", 12);
xlabel("Time", "FontSize", 18);
ylabel("Ratio", "FontSize", 18);

function [t,x] = stoch_sim(x0,mu,beta, t_start, t_end, dt)
    max_iters = (t_end-t_start)/dt;
    
    t = zeros(1,max_iters);
    x = zeros(max_iters,2);
    t(1) = t_start;
    x(1,:) = x0(:);
    for iter = 2:max_iters
        k_i_to_s = mu*dt;
        k_s_to_i = 1- (1-beta*dt).^x(iter-1,2);
        
        x(iter,:) = x(iter-1,:);
        
        n_infected = x(iter,2);
        n_susept = x(iter,1);
        for i = 1:n_infected
            r = rand();
            if r < k_i_to_s
                x(iter,2) = x(iter,2) - 1;
                x(iter,1) = x(iter,1) + 1;
            end
        end
        
        for i = 1:n_susept
            r = rand();
            if r < k_s_to_i
                x(iter,1) = x(iter,1) - 1;
                x(iter,2) = x(iter,2) + 1;
            end
        end
        
        t(iter) = t(iter-1) + dt;
    end
end