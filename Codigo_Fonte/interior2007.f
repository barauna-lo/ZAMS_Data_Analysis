c     interior.for - VERSAO dezembro 2007
c     Curso de Evolucao Estelar 2007B
c     Prof.: Gustavo Porto de Mello
c     Versao revisada do programa de simulacao de interiores estelares.
c     Criado originalmente por Basilio Santiago em 1988
c     Melhorado por diversas geracoes de combatentes...
c
c     Todas as referencias a paginas e tabelas e todas as expressoes
c     utilizadas neste programa estao no livro: Introduction to Stellar
c     Atmospheres and Interiors, E.Novotny 1973, Oxford Univ. Press.
c
c     A partir de agora os comentarios se referem `as linhas de comando
c     do programa, os comentarios sempre precedem a parte do programa a
c     qual se referem.
c
c     Com os comandos abaixo todas as variaveis iniciadas pelas letras
c     a-h,o-z passam a ser consideradas, implicitamente, como variaveis
c     reais com precisao dupla. As variaveis iniciadas por d sao
c     densidades, por vl sao luminosidades, por t temperaturas,
c     por vm massas, por p pressoes, eps se refere `a taxa de prodcao
c     de energia e vkappa `a opacidade.
c
      implicit real*8(a-h,o-z)
      dimension d(400),vl(400),t(400),pg(400),eps(400),vm(400)
      dimension vkappa(400),pr(400), p(400)
c
c     Os parametros a seguir sao dados de entrada que especificam as
c     caracteristicas da estrela que sera' modelada.
c    
      print *, 'Para calcular o modelo de interior estelar entre com'
      print *, 'os seguintes dados da estrela'
      print *, ' ' 
      print *, 'Massa em unidades solares:'
      read(*,*) vm0
      print *, 'Raio em unidades solares:'
      read(*,*) r0
      print *, 'Luminosidade em unidades solares:'
      read(*,*) vl0
      print *, 'Temperatura efetiva em kelvin:'
      read(*,*) t0
      print *, 'Densidade na superficie em g/cm3:'
      read(*,*) d0
      print *, 'Massa percentual de H:'
      read(*,*) x
      print *, 'Massa percentual de metais:'
      read(*,*) z
c     
c     Na primeira parte deste programa definimos todas as constantes que
c     serao usadas nas equacoes do modelo.
c
c     Abaixo sao definidas constantes que serao usadas durante os
c     calculos: rsun, vlsun, vmsun sao raio, luminosidade e massa do Sol
c     respectivamente. r0star, vl0star, vm0star sao raio, luminosidade
c     e massa da estrela que se deseja modelar
c     r0perc sera' usado adiante na definicao do passo da interacao.
c     vk eh a constante de Boltzmann, G eh a constante da gravitacao
c     Atencao, todas as unidades abaixo estao no sistema CGS (g,cm,s)
c
      rsun=6.9598D+10
      vlsun=3.826D+33
      r0star=rsun*r0
      r0perc=0.1D+00*r0star
      vl0star=vlsun*vl0
      vmsun=1.9891D+33
      vm0star=vmsun*vm0
      vk=1.38062D-16
      G=6.67D-08
c
c     Abaixo sao definidas novas constantes, c1 a c9 sao combinacoes das
c     constantes das expressoes para as taxas de geracao de energia
c     pelos ciclos pp e CNO. Essas expressoes estao nas paginas
c     257 e 258 e serao definidas adiante.
c
      pi=3.141593D+00
      ve=2.718281828D+00
      c1=2.36D+10
      c2=-3381.D+00
      c3=1.23D-4
      c4=1.09D-6
      c5=9.5D-10
      c6=7.21D+31
      c7=-1.5231D+04
      c8=3.68D+22
      c9=4.34D+25
      c10=2.52D-15
c
c     Os valores abaixo na realidade nao sao constantes mas para o
c     proposito deste programa serao considerados como constantes.
c     fpp eh o "electron screening factor" ou "electron shielding
c     factor" para o ciclo pp. Este fator depende de t e d e eh dado
c     pela tabela 10-13B (pag. 507). Para o interior do Sol t7 eh 1,6
c     e d = 150, na tabela escolhemos aproximar t7 = 2 e d = 200, pois,
c     aproximar por uma valor maior de t7 deve subestimar fpp e um valor
c     maior de d deve superestima-lo.
c     fcn eh o mesmo para o ciclo CNO. Eh dado pela tabela 10-14B, a
c     mesma aproximacao acima para t7 e d foi usada para manter a
c     consistencia interna.
c     psi eh um fator que leva em conta as diferentes energias efetivas
c     liberadas pelos varios ramos do ciclo pp. A figura 10-7 (pag 507)
c     eh um grafico de psi em funcao de t.
c     gff eh o fator de Gaunt para absorcao livre-livre (pag 131).
c     gbf_t eh a razao entre o fator de Gaunt g_bf e
c     o fator guilhotina t (pag 160) para absorcao ligado-livre (pag 449).
c     Seu valor eh dado em funcao de t, d e x na tabela 10-15 (pag. 509).
c
      fpp=1.04D+00
      fcn=1.33D+00
      psi=1.1D+00
      gff=1.0D+00
      gbf_t=0.31626D+00
c
c     pmm eh o peso molecular medio (adimensional) sob condicao de
c     ionizacao completa.
c     pmmh eh o produto de pmm pela massa do proton
c     a eh a constante de densidade de radiacao, igual a 4*sigma/c
c     (sigma eh a constante de Stephan-Boltzmann).
c     c eh a velocidade da luz, y eh a fracao por massa de He.
c
      x2=x**2.0D+00
      y=1.0D+00-x-z
      pmm=1.0D+00/(2.D+00*x+0.75D+00*y+0.5D+00*z)
      pmmh=1.66053D-24*pmm
      a=7.56471D-15
      c=2.997925D+10
c
c     Definidas as constantes o programa agora executa a primeira
c     iteracao. Primeiramente as variaveis assumem os valores iniciais
c     relativos a camada mais externa do interior estelar, ou seja, elas
c     assumem os valores dos dados de entrada.
c     pg(i) eh a pressao dada pela equacao dos gases perfeitos, pr(i)
c     eh a pressao de radiacao e p(i) eh a pressao total.
c
      d(1)=d0
      vl(1)=vl0star
      vm(1)=vm0star
      t(1)=t0
      pg(1)=d(1)*vk*t(1)/pmmh
      pr(1)=(a*(t0**4))/3
      p(1)=pg(1)+pr(1)
c
c     Agora calcula-se os valores iniciais para as taxas de geracao de
c     energia nuclear pelos ciclos pp (eps1*eps2) e pelo ciclo CNO
c     (eps3) na camada mais externa do interior estelar.
c     As equacoes estao nas paginas 257 e 258. eps eh a taxa total de
c     energia produzida.
c     
      term1=c2*(t(1)**(-1/3.D+00))
      term3=c7*(t(1)**(-1/3.D+00))
      eps1=c1*x2*d(1)*(t(1)**(-2/3.0D+00))*(ve**term1)*psi*fpp
      eps2=1.0D+00+c3*(t(1)**(1/3.0D+00))+c4*(t(1)**(2/3.0D+00))+c5*t(1)
      eps3=c6*x*z*d(1)*(t(1)**(-2/3.0D+00))*(ve**term3)*fcn
      eps(1)=eps1*eps2+eps3
c
c     Calculamos agora a opacidade. vkbf e vkff sao os coeficientes de
c     absorcao medios ligado-livre e livre-livre, respectivamente,
c     sigma_e eh a opacidade devido ao espalhamento eletronico.
c     (ver pag 160) vkappa eh a opacidade total.
c
      vkbf=c9*gbf_t*z*(1.D+00+x)*d(1)*(t(1)**(-3.5D+00))
      vkff=c8*gff*(1.D+00-z)*(1.D+00+x)*d(1)*(t(1)**(-3.5D+00))
      sigma_e=0.2D+00*(1.D+00+x)
      vkappa(1)=vkbf+vkff+sigma_e
c
c     Com todos os valores iniciais (da camada mais externa) definidos o
c     programa ira agora iniciar as proximas iteracoes. Definimos o
c     passo h abaixo, note que o passo eh definido em termos do raio. O
c     programa executa os calculos se aprofundando no raio estelar a
c     partir da camada externa cujos valores foram calculados acima.
c     Abaixo tambem eh aberto o arquivo de saida: "pontos.dat".
c
      i=1
      k=1
      h=r0star/200.D+00
      r=r0star+h/2.D+00
      open(1,file='pontos1.dat')
      open(2,file='pontos2.dat')
      open(3,file='pontos3.dat')
      write(1,11)
 11   format('#      raio    massa/M    densidade    pressao_t'
     &'     Temperatura   Lumin  transporte')
      write(2,21)
 21   format('#  raio       pgas        prad      epspp      epscno'
     &'   epstot')
      write(3,31)
 31   format('#  raio      kapbf      kapff      sigmae      kaptot'
     &'     densidade')
c
c     Abaixo temos uma expressao de relacao, se r for menor que r0perc
c     (10% do raio) o comando sequinte sera executado, caso contrario o
c     comando eh ignorado e o programa segue para a proxima expressao. O
c     indice i diz respeito ao numero da camada que se esta calculando e
c     o indice k eh usado para calcular os gradientes das quantidades de
c     interesse.
c     Calculamos os gradientes radiativo (pag 233), gradt1, e adiabatico
c     (gradt2) de temperatura e escolhemos entre o transporte radiativo
c     ou convectivo. A variavel igrad foi definida para especificar o
c     tipo de transporte de energia em cada camada, aqui 1 significa
c     transporte radiativo e 2 convectivo.
c
10    continue
      if(r.lt.r0perc) h=r0star/400.D+00
      r=r-h
      k=i
      i=i+1
      r2=r**2.D+00
      termo=t(k)**3.D+00
      gradt1=-3.D+00*vkappa(k)*d(k)*vl(k)/(16.D+00*pi*a*c*termo*r2)
      gradt2=-2.D+00*G*d(k)*t(k)*vm(k)/(5.D+00*p(k)*r2)
      gradt1=abs(gradt1)
      gradt2=abs(gradt2)
      if (gradt1.gt.gradt2) go to 20
      gradt=gradt1   
      igrad=1
      go to 30
20    gradt=gradt2
      igrad=2      
c
c     Calculamos agora os incrementos nas variaveis. gradm eh a equacao
c     de continuidade de massa, gradpg o equilibrio hidrostatico, gradpr
c     eh o gradiente da pressao de radiacao, gradl eh o equilibrio
c     termico. Como nosso elemento de raio eh o passo h, estes
c     gradientes devem ser multiplicados por h para obtermos a variacao
c     das grandezas associadas. Somando os gradientes aos valores
c     anteriores obtemos os valores para a proxima camada.
c

30    gradm=4.D+00*pi*r2*d(k)
      gradpg=-G*vm(k)*d(k)/r2
      gradpr=(a*4*(t(k)**3)*gradt)/3
      gradl=4*pi*r2*d(k)*eps(k)
      gradm=abs(gradm)
      gradpg=abs(gradpg)
      gradpr=abs(gradpr)
      gradp=gradpg+gradpr
      gradl=abs(gradl)
      vm(i)=vm(k)-gradm*h
      p(i)=p(k)+gradp*h
      vl(i)=vl(k)-gradl*h
      t(i)=t(k)+gradt*h
      d(i)=(pmmh*p(i))/(vk*t(i))
      pr(i)=c10*(t(i)**4)
c      
      term1=c2*(t(i)**(-1/3.D+00))
      term3=c7*(t(i)**(-1/3.D+00))
      eps1=c1*x2*d(i)*(t(i)**(-2/3.D+00))*(ve**term1)*psi*fpp
      eps2=1.D+00+c3*(t(i)**(1/3.D+00))+c4*(t(i)**(2/3.D+00))+c5*t(i)
      eps3=c6*x*z*d(i)*(t(i)**(-2/3.D+00))*(ve**term3)*fcn
      eps(i)=eps1*eps2+eps3
c 
      vkbf=c9*gbf_t*z*(1.D+00+x)*d(i)*(t(i)**(-3.5D+00))
      vkff=c8*gff*(1.D+00-z)*(1.D+00+x)*d(i)*(t(i)**(-3.5D+00))
      sigma_e=0.2D+00*(1.D+00+x)
      vkappa(i)=vkbf+vkff+sigma_e
c      
c     As expressoes abaixo convertem os parametros r, m e L para
c     unidades percentuais em relacao a estrela. Esses parametros, assim
c     como pressao, densidade e temperatura, darao origem a um arquivo
c     de saida: 'pontos1.dat'.
c
      raio1=r-h/2
      radius=raio1/(rsun*r0)
      vmass=vm(i)/(vmsun*vm0)
      pressg=p(i)
      pressrad= pr(i)
      vlumin=vl(i)/(vlsun*vl0)
      dens=d(i)
      tempt=t(i)
      vkappat=vkappa(i)

c
      write(1,100) radius,vmass,dens,pressg,tempt,vlumin,igrad
100   format(2(1x,f10.3),3(1x,D12.2),1x,f10.3,3x,i3)
200   format(6(1x,D10.3))
300   format(5(1x,D10.3),1(1x,D11.3))

c
      write(2,200)radius,pressg,pressrad,eps1*eps2,eps3,eps1*eps2+eps3
      write(3,300)radius,vkbf,vkff,sigma_e,vkappat,dens
      if(raio1.ge.h) go to 10

      continue
      close (1)
      close(2)
      close(3)
      stop
      end 
