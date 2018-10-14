!########################################################################################################
!	Modulo que faz a dinamica
!########################################################################################################

module dynamics
	use geraRede
	use mod_rndgen
	implicit none
	
	!################################################################################################
	!	Parametros internos ao programa
	!################################################################################################

	integer, private, parameter :: dp = kind(0.0d0)
	
	!################################################################################################
	!	Listas dinamicas
	!################################################################################################
	
	integer, allocatable :: infList(:), susList(:), infSusList(:), sigma(:)


	!################################################################################################
	!	Variaveis dinamicas globais
	!################################################################################################

	
	real(dp), private, parameter :: mu = 1._dp
	
	real(dp) :: lambda
	

	!################################################################################################
	!	Para inicializar
	!################################################################################################
	
	real(dp) :: nInf0
	
	!################################################################################################
	!	Variaveis dinamicas
	!################################################################################################
	
	integer :: nInf, nSus, nInfSus

	
	!################################################################################################
	!	Taxas de eventos
	!################################################################################################

	real(dp) :: rateTotal, lambdaTotal, muTotal

	
	!################################################################################################
	!	Probabilidades associadas aos eventos independentes
	!################################################################################################
	
	real(dp) :: m, l
	

	!################################################################################################
	!	Variaveis temporais
	!################################################################################################

	real(dp) :: t, dt
	
	integer	:: t_pos, tMedio, t_relax, tMax, t_MaxVis
	real(dp) :: t_LS
		
	real(dp), allocatable :: t_medio(:)


	!################################################################################################
	!	Variaveis temporais auxiliares
	!################################################################################################

	
	integer, allocatable :: t_samp(:), samp_surv(:)


	!################################################################################################
	!	Variaveis epidemicas medias
	!################################################################################################

	
	real(dp), allocatable :: rho_medio(:)
	real(dp) :: rho_medioQS, rho2_medioQS
	real(dp) :: Xi


	!################################################################################################
	!	gerador de numeros aleatorios
	!################################################################################################
	
	type(rndgen) :: gen

	
	!################################################################################################
	!	variavel aleatoria
	!################################################################################################
	
	real(dp) :: prob
	
	
	!################################################################################################
	!	Variaveis auxiliares
	!################################################################################################
	integer, allocatable :: listAux(:)
	integer :: nSamples


	!################################################################################################
	!	Vetor que vai guardar a Probabilidade quasi-estacionaria.
	!################################################################################################

	real(dp), allocatable :: Pn_QS(:)


	contains
	
!####################################################################################
!		Aloca listas para calcular medias
!####################################################################################	
	
		subroutine allocaMedias(this, tmax)
			class(grafo), intent(in) :: this
			integer, intent(in) :: tmax
			
			if(allocated(rho_medio)) deallocate(rho_medio)
			allocate(rho_medio(tmax))
			rho_medio = 0._dp

			
			!################################################################################
			!	Lista para guardar medias dos tempos de eventos
			!################################################################################
			
			
!			if(allocated(t_medio)) deallocate(t_medio)
!			allocate(t_medio(tmax))
!			t_medio = 0._dp
			
!			if(allocated(t_samp)) deallocate(t_samp)										! Vai guardar o numero de amostras
!			allocate(t_samp(tmax))															! que participam da dinamica num dado tempo discreto.
!			t_samp = 0
			
!			if(allocated(samp_surv)) deallocate(samp_surv)									
!			allocate(samp_surv(tmax))
!			samp_surv = 0
			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(this%nodes))
				

			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(this%nodes))
						
			t_MaxVis = 0
		end subroutine
			
!####################################################################################
!		Condicao inicial da dinamica
!####################################################################################
		
		
		subroutine condicaoInicial(this)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: i1, j1, k1, lastCand
			
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			class(grafo), intent(in) :: this

			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

			sigma = 0																			! Todos os nos sao suscetiveis
			nInf = int(nInf0 * this%nodes)		 												! Quantidade de nos inicialmente infectados
			nSus = this%nodes - nInf															! Quantidade de nos inicialmente suscetiveis
											
!			write(*,*) "Numero inicial de infectados ", nInf
										
			if(nInf0 < 1._dp)then
				if(allocated(listAux)) deallocate(listAux)	
				allocate(listAux(this%nodes))
			endif
							
			if(nInf0 < 1._dp)then																																				
					do i1 = 1, this%nodes
						listAux(i1) = i1
					enddo
					k1 = 1;	lastCand = this%nodes
					do while(k1 <= nInf)

						i1 = gen%int(1, lastCand)												! lastCand pertence a geraRede
						infList(k1) = listAux(i1)
						sigma(listAux(i1)) = 1

						listAux(i1) = listAux(lastCand)
						lastCand = lastCand - 1

						k1 = k1 + 1
					enddo
			else
					do i1 = 1, this%nodes
						infList(i1) = i1
						sigma(i1) = 1
					enddo
			endif											

			t = 0._dp;
			t_pos = 1;
		end subroutine
		
		
		
		
!####################################################################################	
!		Modelo S-I-S
!####################################################################################
		
		
		subroutine sisProcess(this)


			!###################################################################################
			!	A rede, ou substrato
			!###################################################################################

			class(grafo), intent(in) :: this
			integer ::  n_K
			integer :: i1, j1, k1
			real(dp) :: probInf
			 
			 n_K = 0
			 
			 do i1 = 1, nInf
				n_K = n_K + this%deg(infList(i1))
			 enddo
			 
			 
			
loopDinamico:	do while(t <= tMax)
				
				muTotal = nInf * mu
				lambdaTotal = n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				dt = -1.0_dp * log(1._dp - gen%rnd())/rateTotal
				t = t + dt
				
				prob = gen%rnd()
				if(prob <= m)then																	! Se o processo da vez eh cura.
					i1 = gen%int(1, nInf)
					
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					infList(i1) = infList(nInf)
					nInf = nInf - 1		
!					write(*,*) "Cura! N Infectados = ", nInf
				else	
loop_infecta:				do														! Se o processo da vez eh infeccao.
						i1 = gen%int(1, nInf)
						probInf = 1._dp * this%deg(infList(i1))/this%degMax
						prob = gen%rnd()
						if(prob	<= probInf)then
							j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)
																		! Posicao na lista de adjacencia correspondente aos vizinhos
																											! de um sitio
							if(sigma(this%listAdj(j1)) == 0)then								
								sigma(this%listAdj(j1)) = 1
								n_K = n_K + this%deg(this%listAdj(j1))													! Muda o estado do sitio.
								nInf = nInf + 1													! Numero de infectados aumenta.
								infList(nInf) = this%listAdj(j1)												! Coloco esse novo infectado na primeira posicao disponivel
																											! na lista de infectados.
							endif
		!					write(*,*) "Infeccao! Numero de Infectados= ", nInf
							exit loop_infecta
						endif
					enddo loop_infecta
				endif									
				
				do while(t_pos <= t)
					rho_medio(t_pos) = rho_medio(t_pos) + 1._dp * nInf/this%nodes
					t_medio(t_pos) = t_medio(t_pos) + t
					t_samp(t_pos) = t_samp(t_pos) + 1
					
					if(nInf .ne. 0)then
						t_MaxVis = max(t_MaxVis, t_pos)												! Qual o tempo recorde que uma infeccao
					endif																			! sobreviveu ate agora na rede?	
					t_pos = t_pos + 1
				enddo

				if(nInf == 0) exit loopDinamico														! Porque o estado absorvente foi alcancado

			enddo loopDinamico			
		end subroutine
		
		
!####################################################################################
!		Processo de Contato
!####################################################################################
	
			subroutine contactProcess(this)


			!###################################################################################
			!	A rede, ou substrato
			!###################################################################################

			class(grafo), intent(in) :: this
			integer :: i1, j1, k1
			
			
			m = 1.0_dp * mu/(lambda + mu)
			
loopDinamico:	do while(t <= tMax)
				
				muTotal = nInf * mu
				lambdaTotal = nInf * lambda
				rateTotal = muTotal + lambdaTotal
				
				dt = -1.0_dp * log(1._dp - gen%rnd())/rateTotal
				t = t + dt
				
				prob = gen%rnd()
				if(prob <= m)then																	! Se o processo da vez eh cura.
					i1 = gen%int(1, nInf)
					
					sigma(infList(i1)) = 0
					infList(i1) = infList(nInf)
					nInf = nInf - 1		
!					write(*,*) "Cura! N Infectados = ", nInf
				else																				! Se o processo da vez eh infeccao.
					
					i1 = gen%int(1, nInf)	
					j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)							! Posicao na lista de adjacencia correspondente aos vizinhos
																									! de um sitio
					if(sigma(this%listAdj(j1)) == 0)then								
						sigma(this%listAdj(j1)) = 1													! Muda o estado do sitio.
						nInf = nInf + 1																! Numero de infectados aumenta.
						infList(nInf) = this%listAdj(j1)												! Coloco esse novo infectado na primeira posicao disponivel
																									! na lista de infectados.
					endif
!					write(*,*) "Infeccao! Numero de Infectados= ", nInf
				endif									
				
				do while(t_pos <= t)
					rho_medio(t_pos) = rho_medio(t_pos) + 1.0_dp * nInf/this%nodes
					t_medio(t_pos) = t_medio(t_pos) + t
					t_samp(t_pos) = t_samp(t_pos) + 1
					
					if(nInf .ne. 0)then
						t_MaxVis = max(t_MaxVis, t_pos)												! Qual o tempo recorde que uma infeccao
					endif																			! sobreviveu ate agora na rede?	
					t_pos = t_pos + 1
				enddo

				if(nInf == 0) exit loopDinamico														! Porque o estado absorvente foi alcancado

			enddo loopDinamico			
		end subroutine
						
		
		
!####################################################################################
!			Esta subrotina calcula as medias da epidemia
!####################################################################################
		
		subroutine calculaMedias(this,seed)
			integer :: seed
			integer :: k1, k2
			class(grafo), intent(in) :: this


			!###################################################
			!	Fase de inicializacao da subrotina
			!###################################################
			
			call gen%init(seed)
			
			
			
			do k1 = 1, nSamples
				write(*,*) "Calculando amostra de numero ", k1
				call gastaRnd()
				call condicaoInicial(this)
				call sisProcess(this)
				
			enddo
			
			open(unit=10, file='t_vs_rho.dat', status='unknown')
			do k1 = 1, t_MaxVis
					write(10, *) 1.0_dp * t_medio(k1)/t_samp(k1), ",", 1.0_dp * rho_medio(k1)/nSamples
			enddo
			close(10)

			
			!#############################################
			!	Subrotinas internas
			!#############################################
			
			contains
				subroutine gastaRnd()
					integer :: k4
					real(dp) :: k5
					do k4 = 1, 10 * nSamples
						k5 = gen%rnd()
					enddo
				end subroutine	
				
		end subroutine


!####################################################################################
!		Pede parametros RBS ao usuario
!####################################################################################


		subroutine promptUserRBS()
				write(*,*) " Digite a porcentagem de nos infectados inicialmente"
				write(*,*) ""
				read(*,*) nInf0
						
				write(*,*) " Digite o tempo de relaxacao da epidemia"
				read(*,*) t_relax
				write(*,*) " "

				write(*,*) " Digite o numero de intervalos de tempo para calcular as medias da epidemia"
				write(*,*) " "
				read(*,*) tMedio
				write(*,*) " "
				
				tMax = t_relax + tMedio
				
							
		end subroutine


!#####################################################################################
!	Get Arg RBS
!#####################################################################################

		subroutine getArgRBS()
			integer :: N0
			character(len=20) :: leitura, buffer
			
			N0 = iargc()
			
			if(N0 /= 3) stop " nIn0 real <= 1, t_relax integer, tMedio integer"			

				call getarg(1, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) nInf0 			
				
				
				write(*,*) nInf0
				
				call getarg(2, leitura)
				
				buffer = trim(adjustl(leitura))
				
				read(buffer,*) t_relax 
				
				write(*,*) t_relax
				
				call getarg(3, leitura)
				
				buffer = trim(adjustl(leitura))
				
				read(buffer,*) tMedio 
				
				write(*,*) tMedio
				
				tMax = t_relax + tMedio		
		end subroutine


!####################################################################################
!		Pede parametros ao usuario
!####################################################################################


		subroutine promptUser()
				write(*,*) " Digite a porcentagem de nos infectados inicialmente"
				write(*,*) ""
				read(*,*) nInf0
						
				write(*,*) " Digite a taxa de infeccao da epidemia"
				read(*,*) lambda
				write(*,*) " "

				write(*,*) " Digite o numero maximo de intervalos de tempo para a epidemia"
				write(*,*) " "
				read(*,*) tMax
				write(*,*) " "

				write(*,*) " Qual o numero desejado de amostras?"
				write(*,*) " "
				read(*,*) nSamples
		end subroutine
		
		
		
		
		
		
!####################################################################################
!		Processo SIS com condicao de contorno reflexiva
!####################################################################################

		subroutine sisProcessRBS(this)


			!###################################################################################
			!	A rede, ou substrato
			!###################################################################################

			class(grafo), intent(in) :: this
			integer ::  n_K
			integer :: i1, j1, k1
			real(dp) :: probInf, sumPQS
			real(dp) :: startTotal, finishTotal, timeTotal, startRej, finishRej, timeRej, &
            timeRejAc, propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat


			 n_K = 0
		

			 do i1 = 1, nInf
				n_K = n_K + this%deg(infList(i1))
			 enddo
			 
            
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

			Pn_QS = 0.0_dp
			sumPQS = 0.0_dp

            
loopDinamico:	do while(t <= tMax)
				
				muTotal = nInf * mu
				lambdaTotal = n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				dt = -1.0_dp * log(max(1e-12,gen%rnd()))/rateTotal
				t = t + dt
				
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					sumPQS = sumPQS + dt
				endif
				
				prob = gen%rnd()
		   if(prob <= m)then
																										! Se o processo da vez eh cura.
					if(nInf == 1) cycle loopDinamico
					i1 = gen%int(1, nInf)
					
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					infList(i1) = infList(nInf)
					nInf = nInf - 1

		   else	
loop_infecta:		do		
																										! Se o processo da vez eh infeccao.
						i1 = gen%int(1, nInf)
						probInf = 1._dp * this%deg(infList(i1))/this%degMax
						prob = gen%rnd()

						
						if(prob	> probInf)cycle

							j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																						! Posicao na lista de adjacencia correspondente aos vizinhos
																														! de um sitio
							if(sigma(this%listAdj(j1)) == 1) exit loop_infecta								
								sigma(this%listAdj(j1)) = 1
								n_K = n_K + this%deg(this%listAdj(j1))													! Muda o estado do sitio.
								nInf = nInf + 1																			! Numero de infectados aumenta.
								infList(nInf) = this%listAdj(j1)									

								exit loop_infecta
			enddo loop_infecta
		   endif									
		enddo loopDinamico

			Pn_QS = 1.0_dp * Pn_QS / sumPQS
			
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			
			do i1 = 1, this%nodes
				rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
				rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
			enddo
			rho_medioQS = 1.0_dp * rho_medioQS / this%nodes							! aqui n medio eh normalizado e se torna rho
			
			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((this%nodes)**2.0_dp)						! aqui n^2 medio eh normalizado e se torna rho2
			
			
			Xi = 1.0_dp * this%nodes * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS

                        
		end subroutine

end module
