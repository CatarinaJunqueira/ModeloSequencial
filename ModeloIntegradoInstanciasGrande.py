# -*- coding: utf-8 -*-

from __future__ import print_function
from scipy.io import loadmat
from connect_ifttt import email_alert
#from logsexcel import save_log_excel]
import numpy as np
import cplex
#from compare_models import compare_models

def ModeloIntegradoInstanciasGrande(casename):
    global debug
    debug = False
    print('iniciando o modelo')
    data = loadmat(casename)
    C = data['C'][0][0]
    R = data['R'][0][0]
    Patios = data['Patios'].tolist()
    TT =data['TT']
    phi = data['phi'].tolist()
    Npatios = len(Patios)
    P=Npatios+1 # numero de portos
    flag = 0
    for add in range(Npatios): # corrigindo a quantidade de posicoes vazias nos patios
        if Patios[add][0].shape[0]*Patios[add][0].shape[1] - np.count_nonzero(Patios[add][0]) < Patios[add][0].shape[0] - 1:
            # se o numero de posicoes livres for menor do que o numero de uma coluna - 1
            # adicionar uma linha de zeros no topo do patio
            Patios[add][0] = np.vstack([np.zeros((1,Patios[add][0].shape[1])),Patios[add][0]])
            flag = 1
    if flag == 0:
        print('Nao eh necessario rodar de novo esta instancia')
        return # sair da funcao       
    
    omega=[ [] for i in range(Npatios) ] # omega = conjunto dos indices dos conteineres em cada patio
    for i in range(Npatios):
        Patios[i]=Patios[i][0]
        omega[i]=np.extract(Patios[i]!= 0 , Patios[i]).tolist()
        
    for o in range(Npatios):
        for d in range(P):
            if phi[o][d].shape[1] != 0 :
                phi[o][d] = phi[o][d].tolist()[0]         
            else:
                phi[o][d] = []
    
    
    N=[ 0 for i in range(Npatios) ] # N = quantidade de conteineres em cada patio
    for i in range(Npatios):
        N[i]=np.count_nonzero(Patios[i])
    
    T=N
    
    H=[] # H = numero de linhas de cada patio
    for i in range(Npatios):
        H.append(Patios[i].shape[0])
    
    W= [] # W = numero de colunas de cada patio
    for i in range(Npatios):
        W.append(Patios[i].shape[1])
    print('parametros criados')
    
    model = cplex.Cplex()
    start_time = model.get_time()
    model.objective.set_sense(model.objective.sense.minimize)
    
    #------------------------------------------------------------#
    #--------------------  Variaveis  ---------------------------#
    nvar = 0 
    
    model,xnames,nvar = variavel_x(model,omega,N,H,W,T,nvar)
    print('variavel x criada')
    model,bnames,nvar = variavel_b(model,omega,N,H,W,T,nvar)
    print('variavel b criada')
    model,vnames,nvar = variavel_v(model,omega,N,T,nvar)
    print('variavel v criada')
    model,ynames,nvar = variavel_y(model,omega,N,H,W,T,nvar)
    print('variavel y criada')
    model,znames,nvar = variavel_z(model,omega,N,T,R,C,nvar)
    print('variavel z criada')
    model,qnames,nvar = variavel_q(model,N,R,C,nvar)
    print('variavel q criada')
    model,wnames,nvar = variavel_w(model,N,R,C,nvar)
    print('variavel w criada')
    model,unames,nvar = variavel_u(model,N,R,C,nvar)
    print('variaveis criadas')
    
    #------------------------------------------------------------#
    #-------------------  Restricoes  ---------------------------#
    #Restricoes do Patio:
    model = restricao_P0(model,omega,N,H,W,Patios)
    print('restricao P0 criada')
    model = restricao_P1(model,omega,N,H,W,T)
    print('restricao P1 criada')
    model = restricao_P2(model,omega,N,H,W,T)
    print('restricao P2 criada')
    model = restricao_P3(model,omega,N,H,W,T)
    print('restricao P3 criada')
    model = restricao_P6(model,omega,N,H,W,T)
    print('restricao P6 criada')
    model = restricao_P7(model,omega,N,H,W,T)
    print('restricao P7 criada')
    model = restricao_P8(model,omega,N,H,W,T)
    print('restricao P8 criada')
    model = restricao_P9(model,omega,N,H,W,T)
    print('restricao P9 criada')
    model = restricao_P10(model,omega,N,H,W,T)
    print('restricao 18 criada')
    model = restricao_PA(model,omega,N,H,W,T)
    print('restricao PA criada')
    #Restricoes de Integracao:
    model = restricao_I1(model,omega,N,T)
    print('restricao I1 criada')
    model = restricao_I2(model,omega,N,T)
    print('restricao I2 criada')
    model = restricao_I3(model,omega,N,R,C,T)
    print('restricao I3 criada')
    model = restricao_I4(model,omega,N,R,C,T,P)
    print('restricao I4 criada')
    model = restricao_I5(model,omega,N,R,C,T)
    print('restricao I5 criada')
    model = restricao_I6(model,phi,R,C,T,P,N)
    print('restricao I6 criada')
    model = restricao_I7(model,omega,N,R,C,T)
    print('restricao I7 criada')
    model = restricao_I8(model,omega,N,R,C,T,P)
    print('restricao I8 criada')
    model = restricao_I9(model,omega,N,R,C,P)
    print('restricao I9 criada')
    #Restricoes do Navio:
    model = restricao_N1(model,P,R,C,TT)
    print('restricao N1 criada')
    model = restricao_N2(model,R,C,P)
    print('restricao N2 criada')
    model = restricao_N3(model,N,R,C)  
    print('restricao N3 criada')
    model = restricao_N4(model,P,R,C)
    print('restricao N4 criada')
    model = restricao_I10(model,omega,N,R,C,TT)
    print('restricao I10 criada')
    
    print('restricoes criadas')
    
 #   model.write("InstanciaModeloIntegrado_1.lp")
       
    variaveis = model.variables.get_num()
    print('Numero de Variaveis = ',variaveis)
    
    restricoes = model.linear_constraints.get_num()
    print('Numero de Restricoes = ',restricoes)
    
    z = 'No modelo ha %s variaveis e %s restricoes' %(variaveis,restricoes)
        
    modelotxt = casename + 'ModeloIntegrado.txt'    
    out_file = open(modelotxt,'w+')  
   
  #  model.set_results_stream(out_file)
   
    #model.set_results_stream("Resultados_InstanciaModeloIntegrado_II.txt")
    
#----------------------------------------------------------------------#
#--------------------  Parâmetros do Cplex  ---------------------------#
#----------------------------------------------------------------------#  
      
    model.parameters.timelimit.set(432000) #limite de tempo em segundos (6 dias)
   
   # model.parameters.mip.strategy.search.set(1) #Apply traditional branch and cut strategy; disable dynamic search
   
  #  model.parameters.mip.limits.gomorypass.set(2) #Generate Gomory fractional cuts aggressively
    
  #  model.parameters.emphasis.mip.set(1) # (1) Emphasize feasibility over optimality. (0) Balance optimality and feasibility; default
    
  #  model.parameters.barrier.algorithm.set(0) #0	Default setting, 1	Infeasibility-estimate start, 2	Infeasibility-constant start, 3	Standard barrier
    
  #  model.parameters.mip.strategy.probe.set(2) # Sets the amount of probing on variables to be performed before MIP branching. (2) - Aggressive probing level

    model.parameters.threads.set(20)
    
    model.parameters.mip.display.set(2)
  # parameters.mip.strategy.probe

#----------------------------------------------------------------------#
#----------------------------------------------------------------------#          

    print('resolvendo o modelo')
    start_solve_time = model.get_time()
    model.solve()   
  #  model.solution.pool.write("SolucaoInicial_Instancia_1.sol")
  #  write solution.txt
#    x1=model.solution.progress.get_num_iterations()
#    print('\nNumero de Iteracoes = ',x1)
    
#    x2=int(model.solution.progress.get_num_nodes_processed())
    #x2=model.solution.progress.get_num_nodes_processed()
#    print('\nNumber of nodes processed = ',x2)
     
#    x3=model.solution.progress.get_num_phase_one_iterations()
#    print('\nNumber of iterations to find a feasible solution = ',x3)

    end_time = model.get_time()
    solvetime = end_time - start_solve_time
    total_time = end_time - start_time

    #save log save_log_excel
   # save_log_excel(model,casename,solvetime,total_time,'resultados_integrado.xlsx')
        
    out_file.seek(0)
    out_string =out_file.read()
    out_file.close()
    
    print(out_string) 

    print('Duracao = ',solvetime)
    
    status = model.solution.get_status_string()
    fobj = model.solution.get_objective_value()
    Navio = obterNavio(model,R,C,Npatios)
 #   Sol_Patios = []
 #   for p in range(Npatios):
 #       Pt = obterP(model,H[p],W[p],omega[p],T[p], phi[p])
 #       Sol_Patios.append(Pt)
    
    print('\nSolution status = ',status)
    print('Valor da Funcao Objetivo: ',fobj )
    
    Y = ' Solution status: %s  <br> Valor da Funcao Objetivo do Modelo Integrado: %s  <br> Duracao: %s <br>' %(status,fobj,total_time)
   
    email_alert(Y, z, out_string)
    return model

#------------------------------------------------------------------------------------------------------------------#
#--------------------  Variaveis  ---------------------------#
#------------------------------------------------------------------------------------------------------------------#  
def variavel_x(model,omega,N,H,W,T,nvar):
    xnames = []
    global xind
    xind = dict()

    indx = 0
    for o in range(len(N)):
        for i in range(1,W[o]+1):
            for j in range(1,H[o]+1):
                for k in range(1,W[o]+1):
                    for l in range(1,H[o]+1):
                        for n in omega[o]:
                            for t in range(1,T[o]+1):
                                xnames.append('x_'+str(i)+'_'+str(j)+'_'+str(k)+'_'+str(l)+'_'+str(n)+'_'+str(t))
                                xind[(i,j,k,l,n,t)] = indx
                                indx +=1
    nx = len(xnames)
    
    nvar += nx
    
    lb = [0.0]*nx
    ub = [1.0]*nx
    ctypes =[model.variables.type.binary]*nx
    obj = [1.0]*nx
    model.variables.add(obj=obj, lb=lb, ub=ub, types=ctypes,names=xnames)
    
    return model,xnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_b(model,omega,N,H,W,T,nvar):
    bnames = []
    global bind
    bind = dict()

    indx = nvar
    for o in range(len(N)):
        for i in range(1,W[o]+1):
            for j in range(1,H[o]+1):
                for n in omega[o]:
                    for t in range(1,T[o]+1):
                        bnames.append('b_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
                        bind[(i, j, n, t)] = indx
                        indx += 1
    nb = len(bnames)

    nvar += nb
    
    lb = [0.0]*nb
    ub = [1.0]*nb
    ctypes =[model.variables.type.binary]*nb
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=bnames)
    
    return model,bnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
def variavel_v(model,omega,N,T,nvar):
    vnames = []
    global vind
    vind = dict()

    indx = nvar
    for o in range(len(N)):
        for n in omega[o]:
            for t in range(1,T[o]+1):
                vnames.append('v_'+str(n)+'_'+str(t))
                vind[(n, t)] = indx
                indx += 1
    
    nv = len(vnames)
    
    nvar += nv
    
    lb = [0.0]*nv
    ub = [1.0]*nv
    ctypes =[model.variables.type.binary]*nv
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=vnames)
    
    return model,vnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_y(model,omega,N,H,W,T,nvar):
    ynames = []
    global yind
    yind = dict()

    indx = nvar
    
    for o in range(len(N)):
        for i in range(1,W[o]+1):
            for j in range(1,H[o]+1):
                for n in omega[o]:
                    for t in range(1,T[o]+1):
                        ynames.append('y_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
                        yind[(i,j,n, t)] = indx
                        indx += 1
    ny = len(ynames)
    
    nvar += ny
    
    lb = [0.0]*ny
    ub = [1.0]*ny
    ctypes =[model.variables.type.binary]*ny
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=ynames)
    
    return model,ynames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_z(model,omega,N,T,R,C,nvar):
    znames = []
    global zind
    zind = dict()

    indx = nvar
    
    for o in range(len(N)):
        for n in omega[o]:
            for t in range(1,T[o]+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        znames.append('z_'+str(n)+'_'+str(t)+'_'+str(r)+'_'+str(c))
                        zind[( n, t,r,c)] = indx
                        indx += 1
    nz = len(znames)
    
    nvar += nz
    
    lb = [0.0]*nz
    ub = [1.0]*nz
    ctypes =[model.variables.type.binary]*nz
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=znames)
    
    return model,znames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_q(model,N,R,C,nvar):
    qnames = []
    global qind
    qind = dict()

    indx = nvar
    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
              for r in range(1,R+1):
                  for c in range(1,C+1):
                      qnames.append('q_'+str(o)+'_'+str(d)+'_'+str(r)+'_'+str(c))
                      qind[(o,d, r, c)] = indx
                      indx += 1
    nq = len(qnames)
    
    nvar += nq
    
    lb = [0.0]*nq
    ub = [1.0]*nq
    ctypes =[model.variables.type.binary]*nq
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=qnames)
    
    return model,qnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_w(model,N,R,C,nvar):
    wnames = []                        
    obj = []

    global wind
    wind = dict()

    indx = nvar

    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
            for a in range(o+1,d+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        
                        wnames.append('w_'+str(o)+'_'+str(d)+'_'+str(a)+'_'+str(r)+'_'+str(c))
                        wind[(o, d,a, r, c)] = indx
                        indx += 1
                        if a == d:
                            obj.append(0.0)
                        else:
                            obj.append(1.0)
    
    nw = len(wnames)
    
    nvar += nw
    
    lb = [0.0]*nw
    ub = [1.0]*nw
    ctypes =[model.variables.type.binary]*nw
    model.variables.add(obj=obj, lb=lb, ub=ub, types=ctypes,names=wnames)
    
    return model,wnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_u(model,N,R,C,nvar):
    unames = []
    global uind
    uind = dict()

    indx = nvar
    
    for o in range(1,len(N)+1):  
        for r in range(1,R+1):
            for c in range(1,C+1):
                unames.append('u_'+str(o)+'_'+str(r)+'_'+str(c))
                uind[(o, r, c)] = indx
                indx += 1
    nu = len(unames)
    
    nvar += nu
    
    lb = [0.0]*nu
    ub = [1.0]*nu
    ctypes =[model.variables.type.binary]*nu
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=unames)
    
    return model,unames,nvar


#------------------------------------------------------------------------------------------------------------------#
#-------------------  Restricoes  ---------------------------#
#------------------------------------------------------------------------------------------------------------------#
def restricao_P0(model,omega,N,H,W,Patios):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for n in omega[o]:
                    rest_names.append('restP0_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(o+1))

                    if Patios[o][H[o]-j-1,i] == n :
                        rhs.append(1.0)
                    else :
                        rhs.append(0.0)

                    #add b coeficientes
                    cols = [bind[(i+1,j+1,n,1)]]
                    coefs = [1.0]
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 1: In each time period, each block must either be within the stack or in the outside region:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P1(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(T[o]):
                rest_names.append('restP1_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(1.0)
                cols = [vind[(n,t+1)]]
                coefs = [1.0]
                
                for i in range(W[o]):  
                    for j in range(H[o]):
                        #add b coeficientes
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(1.0)
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 2: In each time period, each slot (i,j) must be occupied by at most one block:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P2(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]):
                for t in range(T[o]):
                    rest_names.append('restP2_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
#Constraint 3: garante que não hajam ‘buracos’ no pátio ao restringir que se há um contêiner posição $(i,j+1)$, 
#então a posição $(i,j)$ abaixo também deve estar ocupada:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P3(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]-1):
                for t in range(T[o]):
                    rest_names.append('restP3_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append(bind[(i+1,j+2,n,t+1)])
                        coefs.append(1.0)
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(-1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 4: restrição de equilíbrio de fluxo entre as variáveis de configuração e de movimento no pátio.  
#Vincula o layout no período t com o layout no período t + 1 através das retiradas e realocações executadas:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P6(model,omega,N,H,W,T):

    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for i in range(W[o]):
            for j in range(H[o]):
                for n in omega[o]:
                    for t in range(1,T[o]):
                        if debug:
                            rest_names.append(
                                'restP6_' + str(i + 1) + '_' + str(j + 1) + '_' + str(n) + '_' + str(t + 1) + '_' + str(o + 1))
                        rhs.append(0.0)
                        cols = [bind[(i + 1, j + 1, n, t + 1)]]  # b_{i,j,n,t}
                        coefs = [1.0]
                        cols.append(bind[(i + 1, j + 1, n, t)])  # -b_{i,j,n,t-1}
                        coefs.append(-1.0)
                        cols.append(yind[(i + 1, j + 1, n, t)])  # y_{i,j,n,t-1}
                        coefs.append(1.0)
                        for k in range(W[o]):
                            for l in range(H[o]):
                                if k != i or l != j:
                                    cols.append(xind[(k + 1, l + 1, i + 1, j + 1, n,t)])  # -\sum_{k=1}^{W}\sum_{l=1}^{H}x_{k,l,i,j,n,t-1}
                                    coefs.append(-1.0)
                                    cols.append(xind[(i + 1, j + 1, k + 1, l + 1, n, t)])  # \sum_{k=1}^{W}\sum_{l=1}^{H}x_{i,j,k,l,n,t-1}
                                    coefs.append(1.0)

                        expr.append(cplex.SparsePair(cols, coefs))


    senses=["E"]*len(expr)

    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)
    return model




#------------------------------------------------------------------------------------------------------------------#
#Constraint 5:  define a variável $v_{nt}$ e assegura que todos os contêineres sejam retirados do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P7(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(1,T[o]):
                rest_names.append('restP7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)
                cols = [vind[(n,t+1)]] #v_{nt}
                coefs = [1.0]                       
                for i in range(W[o]):          
                    for j in range(H[o]):
                        for tt in range(t): #\sum_{k=1}^{W}\sum_{l=1}^{H}x_{k,l,i,j,n,t-1}
                            cols.append(yind[(i+1,j+1,n,tt+1)]) #-\sum_{i=1}^{W}\sum_{j=1}^{H}\sum_{t'=1}^{t-1}y_{ijnt'}
                            coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 6:  garante a política LIFO, ou seja, se no período t, o contêiner $n$ está abaixo do contêiner $q$  
# e o contêiner $n$ é remanejado, então no período $t + 1$ o contêiner $n$ não pode estar alocado em uma posição 
# abaixo do contêiner $q$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P8(model,omega,N,H,W,T):
        
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):

        M = N[o]*((H[o]-1)**2)

        for i in range(W[o]): 
            for k in range(W[o]):
                for j in range(H[o]-1):
                    for l in range(H[o]-1):
                        for t in range(T[o]):
                            rest_names.append('restP8_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(M)
                            cols = [] 
                            coefs = []
                            
                            for n in omega[o]: 
                                cols.append(xind[(i+1,j+1,k+1,l+1,n,t+1)]) #\sum_{n=1}^{N}x_{i,j,k,l,n,t}
                                coefs.append(M)
                            
                                for jj in range(j+1,H[o]):
                                    for ll in range(l+1,H[o]):
                                        cols.append(xind[(i+1,jj+1,k+1,ll+1,n,t+1)])
                                        coefs.append(1.0)               
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 7:  garante que sejam remanejados apenas os contêineres que estão acima, ou seja, na mesma coluna, 
# de um contêiner a ser retirado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P9(model,omega,N,H,W,T):

        
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        M = (H[o] ** 2) * W[o]*(W[o]-1)
        for i in range(W[o]):
            for t in range(T[o]):
                for n in omega[o]: 
                    rest_names.append('restP9_'+str(n)+'_'+str(i+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(M)
                    cols = [] 
                    coefs = []
                    
                    for j in range(H[o]):
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(M)
                        for k in range(W[o]):
                            for l in range(H[o]):
                                for ii in range(i):
                                    cols.append(xind[(ii+1,j+1,k+1,l+1,n,t+1)])
                                    coefs.append(1.0)
                                for iii in range(i+1,W[o]):
                                    cols.append(xind[(iii+1,j+1,k+1,l+1,n,t+1)])
                                    coefs.append(1.0)          
                        
                    expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 8:  garante que nenhum contêiner pode ser remanejado para outra posição que esteja da mesma coluna na
# qual ele se encontra:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P10(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for l in range(H[o]):
                    for n in omega[o]: 
                        for t in range(T[o]):
                            rest_names.append('restP10_'+str(i+1)+'_'+str(j+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(0.0)
                            cols= [xind[(i+1,j+1,i+1,l+1,n,t+1)]]
                            coefs = [1.0]        
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 9:  garante que um contêiner na posição $(i,j)$ só pode ser movido depois que o contêiner na posição 
# $(i,j+1)$ é movido. Se o contêiner na posição $(i,j+1)$ não é movido então temos que $b_{i(j+1)nt} = 1$ e 
# $x_{i(j+1)klnt} = 0$, e o lado esquerdo da equação se torna 0. Consequentemente o lado direito da equação também 
# deve ser 0. Dessa forma, nenhuma realocação ou remanejamento é permitido para o contêiner na posição $(i,j)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_PA(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]-1):
                for t in range(T[o]):
                    rest_names.append('restA_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]:
                        cols.append(bind[(i+1,j+2,n,t+1)]) #\sum_{n=1}^{N}b_{i,j+1,n,t}
                        coefs.append(1.0)
                        cols.append(yind[(i+1,j+1,n,t+1)]) #\sum_{n=1}^{N}y_{i,j,n,t}
                        coefs.append(1.0)
                        for k in range(W[o]):
                            for l in range(H[o]):
                                cols.append(xind[(i+1,j+2,k+1,l+1,n,t+1)]) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j+1,k,l,n,t}
                                coefs.append(-1.0)  
                                cols.append(xind[(i+1,j+1,k+1,l+1,n,t+1)]) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j,k,l,n,t}
                                coefs.append(1.0)

                    expr.append(cplex.SparsePair(cols,coefs))  

    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
# Constraint 10: garante que em cada período de tempo um contêiner seja retirado do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I1(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for t in range(T[o]):
            rest_names.append('restI1_'+str(t+1)+'_'+str(o+1))
            rhs.append(t)
            cols = []
            coefs = []
            for n in omega[o]:          
                #add v coeficientes
                cols.append(vind[(n,t+1)])
                coefs.append(1.0)
                        
            expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model



#------------------------------------------------------------------------------------------------------------------#
# Constraint 11: define a variável $v_{nt}$. Quando um contêiner $n$ é retirado do pátio, a variável $v_{nt}$ se 
# torna 1 e se mantém igual a 1 nos períodos de tempo seguintes:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I2(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                rest_names.append('restI2_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = [vind[(n,t+1)]]
                coefs = [1.0]
                cols.append(vind[(n,t+2)])
                coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model
#------------------------------------------------------------------------------------------------------------------#
# Constraint 12: garante que o contêiner $n$ seja carregado no navio no período de tempo $t$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I3(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                rest_names.append('restI3_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = [vind[(n,t+2)]]
                coefs = [-1.0]
                for r in range(R):
                    for c in range(C):
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(1.0)
                    
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 13: assegura que uma posição $(r,c)$ no navio só pode ser ocupada por um contêiner, seja ele um 
# contêiner que foi carregado no porto atual (porto $o$), em algum porto anterior (porto $o-1$) ou um contêiner
# que já estava no navio e está sendo remanejado em $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I4(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R):
            for c in range(C):
                for t in range(T[o]):
                    rest_names.append('restI4_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append(wind[(oo+1,d+1,a+1,r+1,c+1)])
                                coefs.append(1.0)
                    
                    for d in range(o+1,P):
                        cols.append(qind[(o+1,d+1,r+1,c+1)])
                        coefs.append(1.0)                            
                    
                    for n in omega[o]: 
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(1.0)
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 14: certifica que o contêiner $n$, depois de carregado, não mude de posição enquanto o navio estiver
# parado no mesmo porto:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I5(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                for r in range(R):
                    for c in range(C):
                        rest_names.append('restI5_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1))
                        rhs.append(0.0)     
                       #add z coeficientes
                        cols=[zind[(n,t+1,r+1,c+1)]]
                        coefs = [1.0]
                        cols.append(zind[(n,t+2,r+1,c+1)])
                        coefs.append(-1.0)
                        
                        expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 15:  garante que se há um contêiner na posição $(r,c)$ do navio, ele deve ser um contêiner que acabou 
# de ser embarcado, ou um contêiner de remanejamento:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I6(model,phi,R,C,T,P,N):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for r in range(R):
              for c in range(C):
                  for d in range(o+1,P):
                      rest_names.append('restI6_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1)+'_'+str(d+1))
                      rhs.append(0.0)     
                      #add q coeficientes
                      cols = [qind[(o+1,d+1,r+1,c+1)]]
                      coefs = [1.0]
                      for n in phi[o][d]:
                          cols.append(zind[(n,T[o],r+1,c+1)])
                          coefs.append(1.0)                       
                      for a in range(o+1,d+1):
                          cols.append(wind[(o+1,d+1,a+1,r+1,c+1)])
                          coefs.append(-1.0)                   
                        
                      expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 16: assegura que todos os $N_{o}$ contêineres do pátio $o$ já foram embarcados no navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I7(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for n in omega[o]:
              for t in range(T[o]-1):
                      rest_names.append('restI7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                      rhs.append(1.0)
                      cols = [] 
                      coefs = []
                      for r in range(R):
                          for c in range(C):
                              cols.append(zind[(n,T[o],r+1,c+1)])
                              coefs.append(1.0)                                     
                        
                      expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 17:garante que, durante o processo de carregamento do navio, nenhum contêiner seja alocado em uma 
# posição flutuante ou que ocupe a posição de um contêiner que já estava no navio ou foi remanejado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I8(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                for t in range(T[o]):
                    rest_names.append('restI8_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]: 
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(-1.0)
                        cols.append(zind[(n,t+1,r+2,c+1)])
                        coefs.append(1.0)                      
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append(wind[(oo+1,d+1,a+1,r+1,c+1)])
                                coefs.append(-1.0)
                    
                    for d in range(o+1,P):
                        cols.append(qind[(o+1,d+1,r+1,c+1)])
                        coefs.append(-1.0)                            
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 18: contabiliza o número total de contêineres que foram remanejados no porto $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I9(model,omega,N,R,C,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for d in range(o+1,P):
            rest_names.append('restI9_'+str(o+1)+'_'+str(d+1))
            rhs.append(0.0)
            cols = [] 
            coefs = []
            for oo in range(o):
                for r in range(R):
                    for c in range(C):
                        cols.append(wind[(oo+1,d+1,o+1,r+1,c+1)])
                        coefs.append(1.0)
            for r in range(R):
                for c in range(C):
                    cols.append(qind[(o+1,d+1,r+1,c+1)])
                    coefs.append(-1.0)
                            
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 19: mantém a estabilidade do navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I10(model, omega, N, R, C, TT):
    rhs = []
    rest_names = []
    expr = []

    v = np.sum(TT, axis=0)
    theta = [0.0] * len(N)
    theta[0] = N[0] - v[0]

    for i in range(1, len(N)):
        theta[i] = theta[i - 1] + N[i] - v[i]

    for o in range(len(N)):
        temp = int(np.ceil(float(theta[o]) / float(C)))

        if temp < R:
            rest_names.append('restI10_' + str(o + 1))
            rhs.append(0.0)
            cols = []
            coefs = []

            for c in range(C):
                for r in range(temp, R):
                    cols.append('u_' + str(o + 1) + '_' + str(r + 1) + '_' + str(c + 1))
                    coefs.append(1.0)

            expr.append(cplex.SparsePair(cols, coefs))

    senses = ["E"] * len(expr)

    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs, names=rest_names)
    return model
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 20: restrição de conservação de fluxo, e indica que o número total de contêineres no porto $o$ deve ser
# igual ao número de contêineres que foram embarcados nos portos $p=1,2,...,o$ menos os contêineres que foram 
# desembarcados nos portos $p=1,2,...,o$ :
#------------------------------------------------------------------------------------------------------------------#
def restricao_N1(model,P,R,C,TT):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for d in range(o+1,P):
            rest_names.append('restN1_'+str(o+1)+'_'+str(d+1))
            rhs.append(float(TT[o,d]))
            cols = [] 
            coefs = []
            for a in range(o+1,d+1):        
                for r in range(R):
                    for c in range(C):
                        cols.append(wind[(o+1,d+1,a+1,r+1,c+1)])
                        coefs.append(1.0)
            
            for m in range(o):
                for r in range(R):
                    for c in range(C):            
                        cols.append(wind[(m+1,d+1,o+1,r+1,c+1)])
                        coefs.append(-1.0)                                
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model  
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 21: garante que cada posição $(r, c)$ tenha no máximo um único contêiner:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N2(model,R,C,P):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for r in range(R):
            for c in range(C):
                rest_names.append('restN2_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = [uind[(o+1,r+1,c+1)]]
                coefs = [-1.0]
                for m in range(o+1):
                    for d in range(o+1,P):
                        for a in range(o+1,d+1):
                            cols.append(wind[(m+1,d+1,a+1,r+1,c+1)])
                            coefs.append(1.0)                                                
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 22: garante que existem contêineres embaixo do contêiner que ocupa a célula $(r, c)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N3(model,N,R,C):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                rest_names.append('restN3_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = [uind[(o+1,r+1,c+1)]]
                coefs = [-1.0]
                cols.append(uind[(o+1,r+2,c+1)])
                coefs.append(1.0)
                                              
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 23: é responsável por definir como um contêiner pode ser desembarcado no porto $d$ ao impor que se
# um contêiner que ocupa a posição $(r, c)$, então ele será desembarcado no porto $d$, se não houver um contêiner
# na posição $(r+1, c)$ acima dele:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N4(model,P,R,C):   
    rhs = []
    rest_names = []
    expr = []
    
    for d in range(1,P):
        for r in range(R-1):
            for c in range(C):
                rest_names.append('restN4_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(1.0)
                cols = [] 
                coefs = []
                for o in range(d):
                    for e in range(d,P):
                        cols.append(wind[(o+1,e+1,d+1,r+1,c+1)])
                        coefs.append(1.0)
                    for u in range(d+1,P):
                        for a in range(d+1,u+1):
                            cols.append(wind[(o+1,u+1,a+1,r+2,c+1)])
                            coefs.append(1.0)           
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model   
#------------------------------------------------------------------------------------------------------------------#
# Obter Navio: Retorna como foi a configuração do navio em cada porto 
#------------------------------------------------------------------------------------------------------------------#

def obterNavio(model,R,C,Npatios):

    Navio = []
    for o in range(Npatios):
        Navio.append(np.zeros((R,C)))
        for r in range(R):
            for c in range(C):
                if model.solution.get_values('u_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1)) == 1 :
                    for k in range(o+1):
                        for j in range(o+1,Npatios+1):
                            for v in range(o+1,j+1):
                                if model.solution.get_values('w_'+str(k+1)+'_'+str(j+1)+'_'+str(v+1)+'_'+str(r+1)+'_'+str(c+1)) == 1 :
                                    Navio[o][r,c]=j+1
                                 
    return Navio                                 
#------------------------------------------------------------------------------------------------------------------#
# Obter Patio: Retorna como foi a movimentacao dos patios em cada porto 
#------------------------------------------------------------------------------------------------------------------#
                              
def obterP(model,H,W,omega,T):

    Pt = []                  
 #   for t in range(T):    
 #       Pt.append(np.zeros((H,W)))
 #       for i in range(W):
 #           for j in range(H):
 #               for n in omega:
 #                  if model.solution.get_values('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)) == 1 :                                                     
 #                      Pt[t][H-j-1,i] == n 
 #   return Pt


    for t in range(1,T+1):
        Pt.append(np.zeros((H,W)))
        for i in range(1,W+1):
            for j in range(1,H+1):
                for n in omega:               
                    if model.solution.get_values('b_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t)) == 1 :
                       Pt[t-1][H-j,i-1] = n 
    return Pt
                    
#------------------------------------------------------------------------------------------------------------------#
# FIM
#------------------------------------------------------------------------------------------------------------------#          
                   
if __name__ == "__main__":
    
    X= [8,9,10,11,12]
    for i in X:
        
        name_instance1 = 'InstanciaModeloIntegrado_' + str(i) + '.mat'
        print(name_instance1)
        _ = ModeloIntegradoInstanciasGrande(name_instance1)