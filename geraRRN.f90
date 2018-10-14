!include 'mod_rndgen_multiple.f90'

!##################################################
!	Modulo que gera a rede
!##################################################

module	geraRede
	use mod_rndgen
	implicit none
	
	integer, private, parameter :: dp = kind(0.0d0)
	
	
	type grafo
		integer :: nodes, edge, degMin, degMax, sumDeg
		integer, allocatable:: deg(:), listAdj(:), aux(:), matriz(:,:)
		real(dp) :: degMean, degDev
		contains
			procedure :: iniciaGrafo
	end type
	
	type,extends(grafo) :: grafoRRN
		contains
			procedure ::	iniciaRRN => iniciaGrafoRRN
			procedure ::	liga => ligaRRN
	end type
	

	!###############################################################
	!		Subrotina que inicializa o grafo
	!###############################################################
	
	contains
	
		subroutine iniciaGrafo(this, N)
			integer, intent(in) :: N
			class(grafo), intent(inout) :: this
		
			this%nodes = N;	this%edge = 0
		
			if(allocated(this%deg)) deallocate(this%deg)
				allocate(this%deg(N)); this%deg = 0
			
			if(allocated(this%aux)) deallocate(this%aux)
				allocate(this%aux(N))
				
			
			!########################################################
			! Nao tem como iniciar matriz, aux e nem listAdj
			! ateh que se conheca a distribuicao de graus
			!########################################################
			
			
		end subroutine
		
		
		
		subroutine misturaGrafos(this, p1, k1, k2)
			real(dp), intent(in) :: p1
			real(dp) :: p2
			integer, intent(in) :: k1, k2
			integer :: N1, N2
			integer :: i1, i2
			class(grafo), intent(inout) :: this
			
			if(p1>1) stop "Nao eh possivel misturar grafos com o valor fornecido para p1"
			if(.not.allocated(this%deg)) stop "Inicie primeiro a subrotina iniciaGrafo."
			
			
			N1 = int(1.0_dp * p1 * this%nodes)
			
			do i1 = 1, N1
				this%deg(i1) = k1
			enddo
			
			do i2 = N1+1, this%nodes
				this%deg(i2) = k2	
			enddo
			
			this%aux(1) = 1
			do i1 = 2, this%nodes
				this%aux(i1) = this%aux(i1-1) + this%deg(i1-1)
			enddo
			
			
			this%sumDeg = sum(this%deg)
			
			If(mod(this%sumDeg, 2)/=0) then
				this%deg(this%nodes) = this%deg(this%nodes) + 1
				this%sumDeg = this%sumDeg + 1
			endif
		
		
			if(allocated(this%listAdj)) deallocate(this%listAdj)
			allocate(this%listAdj(this%sumDeg))
			
			this%degMean = 1._dp * sum(this%deg)/this%nodes
			
			this%sumDeg = sum(this%deg)
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
			
		end subroutine
		
		

		!###############################################################
		!		Setter pra inicializar junto ki.
		!		No inicializador padrao nao eh possivel.
		!###############################################################

		subroutine iniciaGrafoRRN(this, degMean)
			class(grafoRRN) :: this
			integer, intent(in) :: degMean
			integer :: i
			
			this%degMean = degMean
			
			this%aux(1) = 1
			
			do i = 2, this%nodes
				this%aux(i) = this%aux(i-1) + this%degMean
			enddo
			
			do i = 1, this%nodes
				this%deg(i) = degMean
			enddo
			
			this%sumDeg = sum(this%deg)
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
			
			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%sumDeg/2, 2))
!				this%matriz = 0
		end subroutine



		!###############################################################
		!		Subrotina que conecta os nos do grafo RRN
		!###############################################################
	
	
		subroutine ligaRRN(this, seed)
			class(grafoRRN) :: this
			integer, allocatable :: lC(:)
			integer :: i, j, k
			integer :: seed
			integer :: cand1, cand2, lastC, nStubs
			integer, allocatable :: degAux(:)
			logical :: floppou
			type(rndgen) :: gen
			
			call gen%init(seed)
			
			this%sumDeg = sum(this%deg)
			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%sumDeg/2, 2))
			
			
			if(allocated(lC)) deallocate(lC)
				allocate(lC(this%sumDeg))
							
			k = 1
			do i = 1, this%nodes
				do j = 1, this%deg(i)
					lC(k) = i
					k = k + 1
				enddo
			enddo

			
			lastC = size(lC)
			nStubs = lastC/2
			write(*,*) "numero de stubs disponiveis", nStubs
			
			if(allocated(degAux)) deallocate(degAux)
			allocate(degAux(this%nodes))
			degAux = 0
			
			
			
do cand1 = 1, this%nodes

Principal:		do while(degAux(cand1) < this%deg(cand1))
						cand2 = gen%int(1, lastC)
						if(degAux(lC(cand2)) == this%deg(lC(cand2)))cycle Principal
						
						if(lC(cand2) == cand1)then
							cycle Principal
						else
							do k = this%aux(cand1), this%aux(cand1) + degAux(cand1) - 1
								if(this%listAdj(k) == lC(cand2)) cycle Principal
							enddo
						endif
												
						this%listAdj(this%aux(cand1) + degAux(cand1)) = lC(cand2)
						this%listAdj(this%aux(lC(cand2)) + degAux(lC(cand2))) = cand1
						degAux(cand1) = degAux(cand1) + 1
						degAux(lC(cand2)) = degAux(lC(cand2)) + 1
						
						this%edge = this%edge + 1
						this%matriz(this%edge, 1) = cand1; this%matriz(this%edge, 2) = lC(cand2)
						lC(cand2) = lC(lastC); lastC = lastC - 1
			enddo Principal
enddo
			write(*,*) "Sobraram ", nStubs - this%edge," stubs"
			
			this%degMax = maxval(this%deg)
			this%degMin = minval(this%deg)
		end subroutine
		
		subroutine comp_Gigante(this)

			!###############################
			!	amostra de rede
			!###############################
			class(grafo) :: this

			!###############################
			!	Variaveis auxiliares
			!###############################
			integer :: k1, k2, k3, k4, k5
			
			integer :: tamUltComp
			
		end subroutine
		

!##################################################################################################
!	Calcula tamanho da componente gigante
!##################################################################################################
		
subroutine sub_classifica_clusters(this)
	class(grafo), intent(in) :: this
	integer, parameter :: dp = kind(0.0d0)
	real(dp) :: clusterQuot
	integer, allocatable :: clustable(:), clusters(:)
	
	integer :: i1, i2, i3, i4
	
	if(allocated(clustable)) deallocate(clustable)
		allocate(clustable(this%nodes))
		clustable = 1													!Todos os nos sao clusterizaveis a priori
		
	if(allocated(clusters)) deallocate(clusters)
		allocate(clusters(this%nodes))
		clusters = 0													!A priori, nao ha clusters
		

	
	end subroutine sub_classifica_clusters
		
		
end module
