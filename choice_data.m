classdef choice_data
    % choice data for one individual with D decisions with potentially
    % differing max M options available.
    %
    % y ... Dx1 index of chosen alt.
    % x ... DxMxK matrix of characteristics.
    % opt ... DxM indicator matrix of options avail. 
    %
    % generates dx ... Dx (M-1) x K chars of differences to chosen alt. 
    %
    % dbauer, 24.9.2015. 
    properties 
        y
        x
        opt
        dx
        M
        K
        D
        reglab
        optlab
    end %properties
    
    methods 
        function obj = choice_data(y,x,opt,ASC,optlab,reglab)
            
            [D,M,K]= size(x);
            K = size(x,3);
            
            obj.M= M;
            obj.D =D;
            obj.K= K;
            if nargin<6
                for k=1:K
                    reglab{k} = sprintf('X%d',k);
                end
            end
            
            if nargin<5
                for m=1:M
                    optlab{m} = sprintf('OPT:%d',m);
                end;
            end;
            if nargin<4
                ASC = 0;
            end;
            if nargin<3 
                opt = ones(D,M);
            end;
            obj.opt = opt; 
            if ASC
                for j=1:M
                    x(:,j,end+1)=1;
                    reglab{end+1} = sprintf('ASC:%d',j);
                end;
                obj.K= K+M;
            end;
            obj.reglab = reglab;
            obj.optlab = optlab;
            dx = NaN(D,M-1,K);
            for d=1:D
                not_chosen = setdiff(1:M,y(d));
                XX = squeeze(x(d,:,:));
                dx(d,:,:) = XX(not_chosen,:)-ones(M-1,1)*XX(y(d),:);
            end;
            obj.y = y;
            obj.x = x;
            obj.dx = dx; 
        end                
    end % methods
end % classdef.

            
    