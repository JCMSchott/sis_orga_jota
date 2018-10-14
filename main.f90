!##################################################
!	Programa principal
!##################################################

program main
	use mod_rndgen
	use geraRede, only: grafoRRN
	use dynamics !, only: contactProcess
	
	implicit none
	
	!##########################################
	!	Parametros internos ao probrama
	!##########################################	
	integer, parameter :: dp = kind(0.0d0)

	!##########################################
	!	Medidor de tempo computacional
	!##########################################	

	real(dp) :: start, finish, start2, finish2, tempSis, p1

	!##########################################
	!	Substrato ou rede e seus
	!	parametros
	!##########################################	


	type(grafoRRN) :: rede
	
	
	integer :: ki, N, k1, k2
	integer :: i7, j7
	integer :: npts
	real(dp) :: dlambda
	!##########################################
	!	Semente para o gerador de numeros
	!	pseudo aleatorios.
	!##########################################	

	integer :: seed
	character(len=255) :: diretorio, diretorioPai, subDiretorio
	
	!##########################################
	!	Inicializacao de arquivos
	!##########################################	



!	diretorioPai = 'Dinamica_Dia_14_Outubro_2018_3'
	
	
	
	
!	subDiretorio = '/Arquivos_Fonte/'
!	diretorio = trim(adjustl(trim(adjustl(diretorioPai))//trim(adjustl(subdiretorio))))
	diretorioPai = ''
	subDiretorio = ''
	diretorio = diretorio//subdiretorio	
	open(unit=89, file = trim(adjustl(diretorio))//'rho_vs_lambda.dat', status = 'unknown')
	open(unit=90, file = trim(adjustl(diretorio))//'Xi_vs_lambda.dat', status = 'unknown')
	open(unit=91, file = trim(adjustl(diretorio))//'n_vs_QS.dat', status = 'unknown')
	open(unit=92, file = trim(adjustl(diretorio))//'tExec_vs_lambda.dat', status = 'unknown')
	open(unit=93, file = trim(adjustl(diretorio))//'conexoes_Rede.csv', status = 'unknown')
	open(unit=150, file = trim(adjustl(diretorio))//'tempoGastoDinamica.dat', status = 'unknown')

!##############################################
!	Inicializacao de parametros importantes
!##############################################

	N = 10**5
	ki = 4
	seed = 999999999
	p1 = 1.0_dp
	k1 = 3
	k2 = 2
	npts = 5000
	dlambda = 0.0020_dp
	lambda = 0.4_dp
!##############################################
!	Inicializacao da rotina
!##############################################

	

	
	
	call getArgRBS()
!	call promptUserRBS()
	call cpu_time(start)
		call rede%iniciaGrafo(N)
		call misturaGrafos(rede, p1, k1, k2)
!		call rede%iniciaRRN(ki)

		write(*,*) " "
		write(*,*) "Rede inicializada"

		
		call rede%liga(seed)

		write(*,*) " "
		write(*,*) "Numero de nos= ", rede%nodes, "grau por no= ", rede%degMean
		write(*,*) " "
				
		 write(93,*) "target", ",", "source"
		do i7 = 1, rede%edge
			write(93,*) rede%matriz(i7, 1),",", rede%matriz(i7,2)
		enddo

		write(*,*) " "
		write(*,*) "Rede gerada com sucesso"
		write(*,*) " "
		
!		seed = 398287654
		
!		call gen%init(seed)
		
		call allocaMedias(rede, tMax)
		
!		call calculaMedias(rede, seed)
		
		do i7 = 1, npts

			lambda = lambda + dlambda

			call condicaoInicial(rede)
			
			if((lambda > 0.520_dp) .and. (lambda < 0.523_dp)) tMedio = 1.0_dp * tMedio/100
			call cpu_time(start2)
				call sisProcessRBS(rede)
			call cpu_time(finish2)
			
				write(89,*) lambda, rho_medioQS
				write(90,*) lambda, Xi
				write(*,*) lambda, Xi				
			

			tempSis = finish2 - start2		
			write(92,*) lambda, tempSis
			write(*,*) lambda, tempSis
			
			if((lambda > 0.520_dp) .and. (lambda < 0.523_dp))then
				do j7 = 1, rede%nodes
					if(Pn_QS(j7) > 0.0_dp)then
						write(91,*) j7, Pn_QS(j7)
						write(*,*) j7, Pn_QS(j7)
					endif
				enddo
			endif
			write(*,*) "amostra de numero ", i7, " concluida"			
		enddo

	write(*,*) "Grau medio da rede", rede%degMean
	call cpu_time(finish)
	write(*,*) "Tempo de execucao= ", finish - start

	close(89)
	close(90)
	close(91)
	close(92)
	close(93)
    close(150)

end program
