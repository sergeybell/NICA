      program testrun

      INTEGER::M
      REAL::free,ftcoef(1000,2),I(1000),X(1000),F(1000)
      OPEN (UNIT=2,file='harmonic.dat',STATUS='unknown')

      M=200
      DO J=1,M,1
      OPEN (UNIT=1,file='sin.dat',STATUS='unknown')
      READ(1,*) I(J),X(J),F(J)
c      WRITE (*,*) I(J),X(J),F(J)
      END DO
      CLOSE(1)

      CALL foft(F,M,ftcoef,free)

      DO J=1,M/2
      WRITE(2,*) ftcoef(J,1), ftcoef(J,2),free
      END DO

      stop
      end

      subroutine foft(dat,nft,ftcoef,free)
c---------------------------------------
c    dat - входной массив данных удвоенной точности
c    nft - число точек в массиве dat
c    ftcoef - массив коэффициентов разложения Фурье размерности (nft/2 * 2) - на выход
c       Один столбец у него для косинусов, другой для синусов
c       Тут для простоты примем, что nft - четное число, тогда длина ftcoef будет ровно nft/2
c    free - свободный член разложения Фурье (умные дяди программеры стандартно запихивают его
c       в выходной массив, но мы так делать не будем, чтобы все было как в формуле) - на выход
c
      integer :: nft,j,n,m,pos
      real :: dat(nft),trig(nft,2),ftcoef(nft/2,2),free,x,twopi

      twopi = 6.28318530717959


c     cчитаем свободный член как среднее арифметическое:
c----------------------------------------------------------
      x = 0.0
      DO j = 1, nft
      x = x + dat(j)
      free = x/nft
      END DO

c   Инициируем trig(nft,2) - рабочий массив косинусов-синусов:
c----------------------------------------------------------

      trig(1,1) = 1.0d0
      trig(1,2) = 0.0d0
      DO j = 2,nft
      trig(j,1) = cos(twopi*(j-1)/nft)
      trig(j,2) = sin(twopi*(j-1)/nft)
      END DO


      DO m = 1,2
      DO j = 1, nft/2-1
      x = 0.0
      DO n = 1, nft
      pos = mod((n-1)*j,nft)+1
      x = x + trig(pos,m)*dat(n)
      END DO
      ftcoef(j,m) = 2*x/nft
      END DO
      END DO
      ftcoef(nft/2,2) = 0.0d0

      x = 0.0d0
      DO n = 1, nft
      pos = mod((n-1)*nft/2,nft)+1
      x = x + trig(pos,1)*dat(n)
      END DO
      ftcoef(nft/2,1) = 2*x/nft

c   Boт поскольку коэффициент ftcoef(nft/2,2) всегда и везде будет 0, ушлые программеры на
c   его место запихивают свободный член разложения. Но мы так не сделаем.

      return
      end
