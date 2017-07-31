
      #######
      function drbasis(nn::Real,qq::Real)

          #inputs
          n = nn
          x = linspace(0,1,n)
          kset = colon(1,n);
          tset = linspace(0,1,n);
          q = qq
          q_global = qq
          n_global = nn

          #constantse
          E = exp(1)
          I = 1im
          Pi = pi

          const0 = 1
          const1 = -sqrt(3)
          const2 = 2*sqrt(3)
          const3 = -sqrt(5)
          const4 = 6*sqrt(5)
          const5 = - 6*sqrt(5)
          const6 = -sqrt(7)
          const7 = 12*sqrt(7)
          const8 = -30*sqrt(7)
          const9 = 20*sqrt(7)
          const10 = 3
          const11 = -60
          const12 = 270
          const13 = -420
          const14 = 210
          const15 = -sqrt(11)
          const16 = 30*sqrt(11)
          const17 = -210*sqrt(11)
          const18 = 560*sqrt(11)
          const19 = 630*sqrt(11)
          const20 = 252*sqrt(11)


          #funcs
         function Power(x,y)
              if(length(x)<2)
                  if(x==E)
                      ee = E;                      # ee = as.brob(E); CORREGIR
                      yRe = real(y)
                      yIm = complex(0,imag(y))
                      ans = (ee^yRe)*(ee^yIm)
      		            #añadimos esta restriccion
      		            if(ee^yRe == Inf)
      			               if(imag(ee^yIm)==0)   # anadimos el o
      				                   ans = complex(ee^yRe)
      			               else
      				                   ans = complex(ee^yRe)
      			               end
      		            end
      		            #  ---  fin
                  else
                    ans = x.^y
                  end
              else
                ans = x.^y
              end
              ans
          end

          function Complex(x,y)
              complex(x,y)
          end

          function turn_Mo(M,Mo)
              function sign_change(x)
                  if(x[1] <= x[2])
                      -1
                  else
                      1
                  end
              end
              D = vcat(mapslices(t->[norm(t,2)],M.+Mo,1),mapslices(t->[norm(t,2)],M.-Mo,1))
              D=transpose(D)

              d = mapslices(t->[sign_change(t)],D,2)
              d=transpose(d)
              cMo = Mo*( eye(length(d)).*d)
              cMo
          end

          function Outer(x,y,fun)
          	mmm = length(x)
          	nnn = length(y)
          	res = eye(mmm, nnn);
      	    ii=1;jj=1
            for ii in 1:mmm
              for jj in 1:nnn
            			res[ii,jj] = fun(x[ii],y[jj])
            	end
      	    end
      	    res
          end

      #############
      #Eigenvalues
      #############

      	  function eigenvaluesn(n,q)
      	    function Q(x,q)
      	        if (q==1)
      	            Q=1
      	        end
      	        if (q==2)
      	            Q=1/3+2*cos(pi.*x) .^ 2/3
      	        end
      	        if (q==3)
      	            Q=2/15+11*cos(pi*x).^2/15+2*cos(pi*x).^4/15
      	        end
      	        if (q==4)
      	            Q=(17+180*cos(pi*x).^2+114*cos(pi*x).^4+4*cos(pi*x).^6)/315
      	        end
      	        if (q==5)
      	            Q=(62 .+ 1072 .* cos(pi .* x).^ 2 .+ 1452 .* cos(pi*x).^4 .+247 .*cos(pi.*x).^6 .+ 2 .*cos(pi.*x).^8)./2835
      	        end
      	        if (q==6)
      	            Q=(1382 + 35396*cos(pi*x).^2 + 83021*cos(pi*x).^4 + 34096*cos(pi*x).^6 + 2026*cos(pi*x).^8 + 4*cos(pi*x).^10)/155295
      	        end
              Q
      	    end

            function sinc(x)
                sin(pi.*x)./(pi.*x)
            end

            s_vec = (pi*(1+1/2*( (q+1)%2 ) +floor((q-1)/2):(n-q)+1/2*( (q+1)%2 )+floor((q-1)/2))).^(2*q)/n
        		j=1:(2*n)
        	  atten = ((sinc(j./(2*n))).^(2*q)./Q(j./(2*n),q))[(q+1):n]
        	  s_vec = s_vec.*atten
      	    s_vec
      	  end

      	 ev = cat(1,zeros(q),eigenvaluesn(n,q))[1:n]   # help

      #############
      #Eigenvectors
      #############

          if(q==1)
              #NULL solution
              nullspace = ones(1:n_global)
              #ODE solution
              kset = (q_global+1):n_global
              function phi1(t,k)
              	sqrt(2)*cos((-1 + k)*Pi*t)
              end
              M = hcat(nullspace,Outer(tset,kset,phi1))./sqrt(n_global)
              Mo =  qr(M)[1]
              Mo = turn_Mo(M,Mo)
      	      Mo
          end

          if(q==2)
              #NULL solution
              const0 = ones(1:n_global)
              const1 = -sqrt(3)
              const2 = 2*sqrt(3)
              nullspace = hcat(const0,
                               const1 + const2.*tset)
              #ODE solution
              kset = (q_global+1):n_global
              function phi2(t,k)
                  ans = complex(
                      			Power(complex(-1),1 + k)/Power(E,(-1.5 + k)*Pi*(1 - t)) +
                      			Power(E,-((-1.5 + k)*Pi*t)) + sqrt(2)*cos(Pi/4 + (-1.5 + k)*Pi*t)
                      		)
      	          real(ans)
              end
              M = hcat(nullspace,Outer(tset,kset,phi2))./sqrt(n_global)
              Mo =qr(M)[1]
              Mo = turn_Mo(M,Mo)
          end

          if(q==3)
              #NULL solution
      	      const0 = ones(1:n_global)
              nullspace = hcat(const0,
                               const1 .+ const2.*tset,
                               const3 .+ const4.*tset.+const5.*tset.^2)
              #ODE solution
              kset = (q_global+1):n_global
              function phi3(t,k)
                  real(complex(
                  	(sqrt(1.5) .- Complex(0,1)./sqrt(2)).*
                  	(Power(complex(-1),1 .+ k)./
                  	Power(E,Power(complex(-1),1/6).*(-2 .+ k).*Pi.*(1 .- t)) .+
                  	Power(E,-(Power(complex(-1),1/6).*(-2 .+ k).*Pi.*t))) .+
                  	(sqrt(1.5) .+ Complex(0,1)./sqrt(2)).*
                  	(Power(complex(-1),1 .+ k)./
                  	Power(E,(Complex(0,-0.5) .+ sqrt(3)./2).*(-2 .+ k).*Pi.*(1 .- t)) .+
                  	Power(E,-((Complex(0,-0.5) .+ sqrt(3)./2).*(-2 .+ k).*Pi.*t))) .-
                  	sqrt(2).*sin((-2 .+ k).*Pi.*t)
                  ))
              end
      	      M = hcat(nullspace,Outer(tset,kset,phi3))./sqrt(n_global)
              Mo =qr(M)[1]
              Mo = turn_Mo(M,Mo)
          end

          if(q==4)
              #NULL solution
      	      const0 = ones(1:nn)
              nullspace = hcat(const0,
                               const1.+const2.*tset,
                               const3.+const4.*tset.+const5.*tset.^2,
                               const6.+const7.*tset.+const8.*tset.^2 .+ const9.*tset.^3)
              #ODE solution
              kset = (q_global+1):nn
              function phi4(t,k)
                  real(complex(
                  		(1.0 .+ sqrt(2)).*(Power(complex(-1),1.0 .+ k)./Power(E,(-2.5 .+ k).*Pi.*(1.0 .- t)) .+
                  		Power(E,-((-2.5 .+ k) .* Pi .* t))) .+
                  		(1.0 ./ sqrt(2) .- Complex(0,1) .* (1.0 .+ 1.0 ./ sqrt(2))) .*
                  		(Power(complex(-1),1.0 .+ k) ./ Power(E,Power(complex(-1),0.25) .* ( -2.5 .+ k).*Pi.*(1.0 .- t)).+
                  		Power(E,-(Power(complex(-1),0.25) .* (-2.5 .+ k).*Pi.*t))) .+
                  		(1.0 ./ sqrt(2) .+ Complex(0,1) .* (1.0 .+ 1.0 ./ sqrt(2))).*
                  		(Power(complex(-1),1.0 .+ k) ./
                  		Power(E,(Complex(1,-1) .* (-2.5 .+ k) .* Pi .* (1.0 .- t)) ./ sqrt(2)) .+
                  		Power(E,(Complex(-1,1) .* (-2.5 .+ k) .* Pi .* t) ./ sqrt(2))) .-
                  		sqrt(2) .* cos(Pi ./ 4.0 .- (-2.5 .+ k) .* Pi .* t)
                  	))
              end
      	      M = hcat(nullspace,Outer(tset,kset,phi4))./sqrt(n_global)
              Mo =qr(M)[1]
              Mo = turn_Mo(M,Mo)
          end

          if(q==5)
              #NULL solution
      	       const0 = ones(1:n_global)
               nullspace = hcat(const0,
                                const1 .+ const2 .* tset,
      			                    const3 .+ const4 .* tset .+ const5 .* tset .^ 2,
      			                    const6 .+ const7 .* tset .+ const8 .* tset .^ 2 .+ const9 .* tset .^ 3,
      			                    const10 .+ const11 .* tset .+ const12 .* tset .^ 2 .+ const13 .* tset .^ 3 .+ const14 .* tset .^ 4)
              #ODE solution
              kset = (q_global+1):n_global
       	      function phi5(t,k)
                  real(complex(
                  		(sqrt(2) .*(1 .+ sqrt(5) ./ 2.0) .- (Complex(0,0.5) .*
                  		(sqrt(10 .- 2 .* sqrt(5)) .+ sqrt(2 .* (5 .+ sqrt(5))))) ./ sqrt(2)) .*
                  		(Power(complex(-1),1 .+ k) ./ Power(E,Power(complex(-1),0.1) .* (-3 .+ k) .* Pi .* (1 .- t)) .+
                  		Power(E,-(Power(complex(-1),0.1) .* (-3 .+ k) .* Pi .* t))) .+
                  		(-(1 ./ sqrt(2)) .- (Complex(0,0.5) .*
                  		(sqrt(10 .- 2 .* sqrt(5)) .+ sqrt(2 .* (5 .+ sqrt(5))))) ./ sqrt(2)).*
                  		(Power(complex(-1),1 .+ k) ./ Power(E,Power(complex(-1),0.3) .* (-3 .+ k) .* Pi .* (1 .- t)) .+
                  		Power(E,-(Power(complex(-1),0.3) .* (-3 .+ k) .* Pi .* t))) .+
                  		(sqrt(2) .* (1 .+ sqrt(5) ./ 2.0) .+
                  		(Complex(0,0.5) .* (sqrt(10 .- 2 .* sqrt(5)) .+ sqrt(2 .* (5 .+ sqrt(5))))) ./ sqrt(2)) .*
                  		(Power(complex(-1),1 .+ k) ./
                  		Power(E,(sqrt(0.625 .+ sqrt(5) ./ 8.0) .- Complex(0,0.25) .* (-1 .+ sqrt(5))) .* (-3 .+ k) .* Pi .* (1 .- t)) .+
                  		Power(E,-((sqrt(0.625 .+ sqrt(5) ./ 8.0) .- Complex(0,0.25) .* (-1 .+ sqrt(5))) .* (-3 .+ k) .* Pi .* t))) .+
                  		(-(1.0 ./ sqrt(2)) .+ (Complex(0,0.5) .* (sqrt(10 .- 2 .*sqrt(5)) .+ sqrt( 2 .* (5 .+ sqrt(5))))) ./ sqrt(2)).*
                  		(Power(complex(-1),1.0 .+ k) ./
                  		Power(E,(sqrt(0.625 .- sqrt(5) ./ 8.0) .- Complex(0,0.25) .* (1.0 .+ sqrt(5))) .* (-3 .+ k) .* Pi .* (1.0 .- t)) .+
                  		Power(E,-((sqrt(0.625 .- sqrt(5) ./ 8.0) .- Complex(0,0.25) .* (1.0 .+ sqrt(5))).*
                  		(-3.0 .+ k) .* Pi .* t))) .- sqrt(2) .* cos((-3.0 .+ k) .* Pi .*t)
                  		))
              end
      	      M = hcat(nullspace,Outer(tset,kset,phi5))./sqrt(n_global)
              Mo =qr(M)[1]
              Mo = turn_Mo(M,Mo)
          end

          if(q==6)
              #NULL solution
      	       const0 = ones(1:n_global)
               nullspace = hcat(const0,
                                const1 .+ const2 .* tset,
      			                    const3 .+ const4 .* tset .+ const5 .* tset .^ 2,
      			                    const6 .+ const7 .* tset .+ const8 .* tset .^ 2 .+ const9 .* tset .^ 3,
      			                    const10 .+ const11 .* tset .+ const12 .* tset .^ 2 .+ const13 .* tset .^ 3 .+ const14 .* tset .^ 4,
                                const15 .+ const16 .* tset .+ const17 .* tset .^ 2 .+ const18 .* tset .^ 3 .+ const19 .* tset .^ 4 .+ const20 .* tset .^ 5)
              #ODE solution
              kset = (q_global+1):n_global
       	      function phi6(t,k)
                  real(complex(
                  		(3 .+ 2.0 .* sqrt(3)) .*
                      (Power(complex(-1),1 .+ k) ./ Power(E,(-3.5 .+ k) .* Pi .* (1 .- t)) .+
                      Power(E,-((-3.5 .+ k).* Pi .* t))) .+
                      ((1 .+ sqrt(3)) ./ 2.0 .- Complex(0,0.5) .* (5 .+ 3 .* sqrt(3))) .*
                      (Power(complex(-1),1 .+ k) ./ Power(E,Power(complex(-1),1.0 ./ 6) .* (-3.5 .+ k) .* Pi .* (1 .- t)) .+
                      Power(E,-(Power(complex(-1),1.0 ./ 6) .* (-3.5 .+ k) .* Pi .* t))) .+
                      ((-3 .- sqrt(3))./2.0 .- Complex(0,0.5) .* (1 .+ sqrt(3))) .*
                      (Power(complex(-1),1 .+ k) ./
                      Power(E,Power(complex(-1),1.0 ./ 3) .* (-3.5 .+ k) .* Pi .* (1 .- t)) .+
                      Power(E,-(Power(complex(-1),1 ./ 3) .* (-3.5 .+ k) .* Pi .* t))) .+
                      ((-3 .- sqrt(3)) ./ 2.0 .+ Complex(0,0.5) .* (1 .+ sqrt(3))) .*
                      (Power(complex(-1),1 .+ k) ./
                      Power(E,(0.5 .- Complex(0,0.5) .* sqrt(3)) .* (-3.5 .+ k) .* Pi .* (1 .- t)) .+
                      Power(E,-((0.5 .- Complex(0,0.5) .* sqrt(3)) .* (-3.5 .+ k) .* Pi .* t))) .+
                      ((1 .+ sqrt(3)) ./ 2.0 .+ Complex(0,0.5) .* (5 .+ 3 .* sqrt(3))) .*
                      (Power(complex(-1),1 .+ k) ./
                      Power(E,(Complex(0,-0.5) .+ sqrt(3) ./ 2.0) .* (-3.5 .+ k) .* Pi .* (1 .- t)) .+
                      Power(E,-((Complex(0,-0.5) .+ sqrt(3) ./ 2.0) .* (-3.5 .+ k) .* Pi .* t))) .-
                      sqrt(2) .* cos(Pi ./ 4.0 .+ (-3.5 .+ k) .* Pi .* t)
                  		))
              end
      	      M = hcat(nullspace,Outer(tset,kset,phi6))./sqrt(n_global)
              Mo =qr(M)[1]
              Mo = turn_Mo(M,Mo)
          end

        	eigenvalues = ev
        	eigenvectorsQR = Mo
        	eigenvectors = M

      	#final output
      	list = (eigenvectors, eigenvectorsQR, eigenvalues, x)

      end

#####################################################################################################################################################

      function holajulia()
        msm = "hola mundo desde julia-pkg"
        msm
      end


#####################################################################################################################################################
holajulia()
