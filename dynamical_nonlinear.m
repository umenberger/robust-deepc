classdef dynamical_nonlinear
   
    properties (Access = private)
        
        T = 500; % initial internal memory buffer
                      
    end
    
    properties 
       
% parameters
        A = nan;
        B = nan;
        Bw = nan;
        C = nan;
        D = nan;
        Dw = nan;
        
        param_nl_x_1 = nan;
        param_nl_x_3 = nan;
        param_nl_y_1 = nan;
              
        nx = nan;
        nu = nan;
        nw = nan;
        ny = nan;
        
        x1 = nan;
        
        t = 1; % internal time
        
        x = nan; 
        u = nan;
        w = nan;
        y = nan;
                  
    end
    
    
    
    methods
        
%% method: constructor
% for now we just pass in the initial condition 

        function obj = dynamical_nonlinear(x1)
            
            obj.x1 = x1;            
            
            T = obj.T;
            
            nx = 4;
            nu = 1;
            ny = 1;
            nw = 1;
            
            obj.nu = 1;
            obj.nw = 1;
            obj.ny = 1;             
            
            obj.nx = nx;
            obj.nu = nu;
            obj.nw = nw;
            obj.ny = ny;            
            
            obj.x = nan(nx,T); 
            obj.u = nan(nu,T);
            obj.w = nan(nw,T);
            obj.y = nan(ny,T);  
            
            obj.x(:,1) = x1;
            
% default values of parameters            
            A = [0.9, -0.2, 0.8, -0.3;
                 -0.3, -0.8, 0.5, 0.4;
                 0.5, 0.4, 0.7, -0.3;
                 0.2, -0.3, 0.5, -0.8];
            A = real(0.9*A/spectralRadius(A));     
            B = [1; zeros(3,1)];
            Bw = [0.2;0;0;0]; % how disturbances enter   
            
            C = 1;
            D = 0;
            Dw = 0.05;
            
            obj.A = A;
            obj.B = B;
            obj.Bw = Bw;
            obj.C = C;
            obj.D = D;
            obj.Dw = Dw;
                        
            obj.param_nl_x_1 = -0.05;
            obj.param_nl_x_3 = 0;
            obj.param_nl_y_1 = 0.05;            
            
        end
        
%% method: state_transition         
        function x_ = state_transition(obj,x,u,w)           
            
            x_ = obj.A*x + obj.B*u + obj.Bw*w;
            
            x_(1) = x_(1) + obj.param_nl_x_1*obj.nonlinear_spring(x(1));
            x_(3) = x_(3) + obj.param_nl_x_3*obj.nonlinear_spring(x(3) - x(1));
        
        end

        function y = output_map(obj,x,u,w)
            
            y = obj.C*x(1) + 0*x(1).^2 + obj.param_nl_y_1*x(1).^3 + obj.Dw*w;
            
        end
        
        function obj = simulate(obj,u,w)
            
            N = size(u,2);
            
            for i = 1:N
               
                t = obj.t;
                
                ut = u(:,i); % this is kind of confusing...
                wt = w(:,i);
                xt = obj.x(:,t);
                
                obj.u(:,t) = ut;
                obj.w(:,t) = wt;
                
                obj.x(:,t+1) = state_transition(obj,xt,ut,wt);
                
                obj.y(:,t) = output_map(obj,xt,ut,wt);
                
                obj.t = t + 1;
                
            end
            
        end
        
%% get history
        function [x,u,w,y] = get_history(obj)
            
            x = obj.x(:,1:obj.t-1);
            u = obj.u(:,1:obj.t-1);
            w = obj.w(:,1:obj.t-1);
            y = obj.y(:,1:obj.t-1);

        end                  
        
        
    end
    
    
    methods (Static)
    
        function f = nonlinear_spring(s)
            
            f = 0.2*s + 0.1*s.^2 + 0.6*s.^3;
            
        end
               
    end
    
end















