function dim = get_perm_handle(Y,X,data_parameters,simu_parameters)
alpha = simu_parameters.alpha;
h = simu_parameters.nslices;
nsample = simu_parameters.nsamples;
[n,p] = size(X);

dim = @perm;

    function d = perm(Fo,dFo)
        fpo = Fo(eye(p));
        for u=1:p-1,
            guess = get_initial_estimate(Y,X,u,data_parameters,simu_parameters);
            Wo = guess(Fo);
            if ~isempty(simu_parameters.sg),
                [fno,Wu] = sg_min(Fo,dFo,Wo,simu_parameters.sg{:},simu_parameters.maxiter);
            else
                [fno,Wu] = sg_min(Fo,dFo,Wo,'prcg','euclidean',{1:u},'quiet',simu_parameters.maxiter);
            end
            T = 2*(fno-fpo);
 
            [Q R] = qr(Wu);
            Wu0 = Q(:,(u+1):p);
 

            Xnew = X*[Wu Wu0];
            Xrsp = zeros(n,p);
            fns = 1:nsample;
            fps = 1:nsample;
            for i=1:nsample,
                YY = randperm(n); % ver implementación alternativa en denboot.m
                if u > 0,
                    Xrsp(:,1:u) = Xnew(:,1:u);
                end
                Xrsp(:,(u+1):p) = Xnew(YY(:),(u+1):p);
                newpars = setdatapars_v2(Y,Xrsp(:,:),h);
                Fnew = F(@F4lad,newpars);
                dFnew = dF(@dF4lad,newpars);

                fps(i) = Fnew(eye(p));
                    guess = get_initial_estimate(Y,Xrsp(:,:),u,newpars,simu_parameters);
                    Wo = guess(Fnew);
                    if ~isempty(simu_parameters.sg),
                        fns(i) = sg_min(Fnew,dFnew,Wo,simu_parameters.sg{:},simu_parameters.maxiter);
                    else
                        fns(i) = sg_min(Fnew,dFnew,Wo,'prcg','euclidean',{1:u},'quiet',simu_parameters.maxiter);
                    end
            end
            fnmfp = 2*(fns - fps);
            if (sum(fnmfp>T)/nsample > (alpha)) || (u==p-1),
                d = u;
                break;
            end
        end
   end
end

    