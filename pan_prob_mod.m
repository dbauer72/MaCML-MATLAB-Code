classdef pan_prob_mod
    % panel probit model specification
    % of the form 
    % U_{aij} = sum X_{aijk} (beta_k + z_{jk}) + e_{aij} 
    % for alternative a in decision i for person j.
    % The individual specific slopes z_{j:} are assumed drawn from
    % a normal density with zero mean and variance Omega = O*O'.
    % The vector e_{:ij} is N(0,Sigma), Sigma = L*L'.
    %
    % Parameterisation: 
    %     beta = H_beta * theta_beta + f_beta. 
    %     vec(O) = H_O *theta_O + f_O
    %     vec(L)= H_L theta_L + f_L
    %
    % dbauer, 23.9.2015. 
    properties 
        b
        nvarb
        thb
        Hb
        fb
        db
        O
        nvarO
        thO
        HO
        fO
        L
        nvarL
        thL
        HL
        fL
        Sigma
        Omega
        MSigma
        dMSigma
        dOmega
        m
        k 
        Opt_Names
        Reg_Names
    end
    
    methods 
        function obj = pan_prob_mod(thb,Hb,fb,thO,HO,fO,thL,HL,fL,Opt_Names,Reg_Names)
            
            if nargin==2
                theta = thb;
                ppm = Hb; 
                 try 
                     obj = pan_prob_mod(theta(1:ppm.nvarb),ppm.Hb,ppm.fb,theta(ppm.nvarb+[1:ppm.nvarO]),ppm.HO,ppm.fO,theta(ppm.nvarb+ppm.nvarO+[1:ppm.nvarL]),ppm.HL,ppm.fL,ppm.Opt_Names,ppm.Reg_Names);
                 return
                 catch
                     error('Wrong constructor usage!');
                 end
            end
            k = size(Hb,1);
            obj.k = k;
            km = sqrt(size(HO,1));          
            m = sqrt(size(HL,1));   
            obj.m = m;
            if km ~= k
                error('Dimensions of Hb and HO do not match!');           
            end;
            
            %beta
            obj.b=Hb*thb(:)+fb;
            obj.Hb = Hb;
            obj.fb=fb;
            obj.thb = thb;
            obj.nvarb = length(thb);
            obj.db = Hb;
            
            % Omega.
            obj.HO = HO;
            obj.fO = fO;
            obj.thO = thO;
            obj.nvarO = length(thO);
            obj.O = reshape(HO*thO(:)+fO,k,k);
            obj.Omega = obj.O*obj.O';
            obj.dOmega = zeros(k,k,obj.nvarO);
            for j=1:obj.nvarO
                dO = reshape(HO(:,j),k,k);
                obj.dOmega(:,:,j)=  obj.O*dO'+dO* obj.O';
            end;
            
            % Sigma.
            try 
                obj.L = reshape(HL*thL(:)+fL,m,m);
            catch
                error('Dimensions of HL and m do not match!');
            end;
            obj.thL = thL;
            obj.nvarL = length(thL);
            obj.HL =HL;
            obj.fL = fL;
            obj.Sigma = obj.L*obj.L';
            obj.MSigma = zeros(m-1,m-1,m);
            obj.dMSigma = zeros(m-1,m-1,m,length(thL));
            for j=1:m
                    in = setdiff(1:m,j);
                    obj.MSigma(:,:,j)=obj.Sigma(in,in) - ones(m-1,1)*obj.Sigma(j,in)  - obj.Sigma(in,j)*ones(1,m-1)+obj.Sigma(j,j);              
            end
            for l=1:length(thL)
                dL = reshape(HL(:,l),m,m);
                dSigma = dL* obj.L'+ obj.L*dL';
                
                for j=1:m
                    in = setdiff(1:m,j);
                    obj.dMSigma(:,:,j,l)=dSigma(in,in) - ones(m-1,1)*dSigma(j,in)  - dSigma(in,j)*ones(1,m-1)+dSigma(j,j);              
                end
            end
        
            % Names 
            if length(Opt_Names) == m
                obj.Opt_Names = Opt_Names;
            else 
                obj.Opt_Names =cell(0);
                for jm=1:m
                    obj.Opt_Names{end+1}=sprintf('Alt.%d',jm);
                end;
            end;
            if length(Reg_Names) == k
                obj.Reg_Names = Reg_Names;
            else 
                obj.Reg_Names =cell(0);
                for jk=1:k
                    obj.Reg_Names{end+1}=sprintf('X.%d',jk);
                end;
            end;   
        end % end constructor 
        
        function theta = extract_par(mod)
            % returns the vector of parameters.
            theta = [mod.thb(:)',mod.thO(:)',mod.thL(:)'];
        end    
            
    end % end methods
end % end classdef
