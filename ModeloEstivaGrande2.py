# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 12:18:55 2019

@author: catar
"""

from __future__ import print_function
from scipy.io import loadmat
import numpy as np
import cplex
from connect_ifttt import email_alert
import socket
from ModeloSequencial import ModeloSequencial

def ModeloEstivaGrande2(casename):
    
    global debug
    debug = False
    data = loadmat(casename)
    C = data['C'][0][0]
    R = data['R'][0][0]
    Patios = data['Patios'].tolist()
    TT =data['TT']
    Npatios = len(Patios)
    P=Npatios+1 # numero de portos
    for add in range(Npatios): # corrigindo a quantidade de posicoes vazias nos patios
        if Patios[add][0].shape[0]*Patios[add][0].shape[1] - np.count_nonzero(Patios[add][0]) < Patios[add][0].shape[0] - 1:
            # se o numero de posicoes livres for menor do que o numero de uma coluna - 1
            # adicionar uma linha de zeros no topo do patio
            Patios[add][0] = np.vstack([np.zeros((1,Patios[add][0].shape[1])),Patios[add][0]])
    phi = data['phi'].tolist()
    
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
    #------------------------------------------------------------#
    nvar = 0 
    model,wnames,nvar = variavel_w(model,N,R,C,nvar)
    print('variavel w criada')
    model,unames,nvar = variavel_u(model,N,R,C,nvar)
    print('variavel u criada')

    print('variaveis criadas')
    
    #------------------------------------------------------------#
    #-------------------  Restricoes  ---------------------------#    
    model = restricao_I10(model,omega,N,R,C,TT)
    print('restricao I10 criada')    
    model = restricao_N1(model,P,R,C,TT)
    print('restricao N1 criada')    
    model = restricao_N2(model,R,C,P)
    print('restricao N2 criada')
    model = restricao_N3(model,N,R,C)
    print('restricao N3 criada')    
    model = restricao_N4(model,P,R,C)
    print('restricao N4 criada')
    
    print('restricoes criadas')
    
    #model.write(casename+"Estiva.lp")
    
    variaveis = model.variables.get_num()
    print('Numero de Variaveis = ',variaveis)    
    restricoes = model.linear_constraints.get_num()
    print('Numero de Restricoes = ',restricoes)
    
    z = 'No modelo do plano de estiva ha %s variaveis e %s restricoes' %(variaveis,restricoes)
    
    modelotxt = casename + 'Estiva.txt'    
    out_file = open(modelotxt,'w+')  
    model.set_results_stream(out_file)    
    model.set_results_stream(modelotxt)
    
    model.parameters.timelimit.set(43200)  # limite de tempo em segundos (12 horas)
    model.parameters.threads.set(20)
    #model.parameters.output.writelevel.set(1) # para escrever a solução 

    print('resolvendo o modelo')
    model.solve()
    #model.solution.write("lpex.sol")
    
    out_file.seek(0)
    out_string =out_file.read()
    out_file.close()
    
    print(out_string)

    end_time = model.get_time()
    solvetime = end_time-start_time
    print('Duracao = ',solvetime)
    
    status = model.solution.get_status_string()
    fobj = model.solution.get_objective_value()
    #Navio = obterNavio(model,R,C,Npatios)
  
    Y = 'Instancia Estiva: %s <br> Rodado em: %s <br> Solution status: %s  <br> Valor da Funcao Objetivo Carregamento: %s  <br> Duracao: %s <br>' %(casename,socket.gethostname(),status,fobj,solvetime)
    email_alert(Y,z, out_string)
    
    print('\nSolution status = ',status)
    print('Valor da Funcao Objetivo Estiva: ',fobj )
    
    SolucaoCarregamento, fobjCarregamento = ModeloSequencial(casename,model)
    
    fobj2 = fobj + fobjCarregamento
    print('Valor da Funcao Objetivo Sequencial: ',fobj2)
               
#------------------------------------------------------------#
#--------------------  Variaveis  ---------------------------#
#------------------------------------------------------------#

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

#------------------------------------------------------------#
#------------------------------------------------------------#
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

#------------------------------------------------------------#
#-------------------  Restricoes  ---------------------------#
#------------------------------------------------------------#


#------------------------------------------------------------#
# Constraint 19: mantém a estabilidade do navio:
#------------------------------------------------------------#
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
                              
def obterP(model,H,W,omega,T,phi):

    Pt = []                  
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
    
    X = ['7_1', '8_1']
    for i in X:
        
        name_instance1 = 'InstanciaModeloIntegrado_' + i + '.mat'
        print(name_instance1)
        _ = ModeloEstivaGrande2(name_instance1)          