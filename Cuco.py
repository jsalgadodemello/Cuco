# This Python file uses the following encoding: utf-8

import random
import numpy as np
from math import pi, sin, cos, gamma
import os
import time
# coding=<encoding name>
print('\nBem Vindo ao programa de cálculo de configuração para recarga do reator PWR\n')
n = int(input('Número de Ninhos: '))     # Número de Ninhos
max_geracoes = int(input('Número de Gerações: '))   # Número Máximo de gerações que a ser executada pelo algoritmo

t_inicio = time.ctime()
dimensao = 20   # Dimensão do problema (bidimensional para f(x) = somatório de  xi²)
T = 1e-7   # Critério de Parada/ Precisão (10^-7 é próximo de zero)
pa = 0.25   # Probabilidade de Abandono do ninho

# Definindo o espaço de busca
li = -1*np.ones((n, dimensao)) #Limite inferior do espaço de busca
ls = 1*np.ones((n, dimensao)) #Limite superior do espaço de busca

# Criação dos Ninhos
ninho = np.zeros((n, dimensao))
for i in range(n):
    for j in range(dimensao):
        ninho[i][j] = li[i][j] + (ls[i][j] - li[i][j]) * np.random.rand()
    # Substituindo cada posição do ninho por um valor aleatório dentro do espaço de busca
# Definindo valor de Fitness MUITO ruim, para início do cálculo  
# (todos os valores que forem achados depois desse serão melhores)
fitness = 1.e10 * np.ones((n,1))

###################################################################################################
# Lendo arquivo de Saída do Recnod (salvando as variáveis para fator de pico e concentração de boro)
def lendo_boro_pico():
    aux = []
    arquivo = open('Jessica_Pico_Boro.dat', 'r')
    content = arquivo.read().split("\n")
    for i in range(len(content) - 1):
        linha2 = content[i].replace(',', '.')
        aux.append(float(linha2))
    if aux[0] == 10:
        pico = aux[0]
    else:
        pico = aux[0]/1000
    boro = aux[1]
    arquivo.close()

    return pico, boro
###################################################################################################

###################################################################################################
# Criando arquivo de entrada do recnod
def criando_entrada(ninho_x):
    ninho_1 = []
    ninho_2 = []
    posicoes_1 = []
    posicoes_2 = []
    posicoes = np.zeros((n, dimensao))
    for i in range(int(len(ninho_x))):
        ninho_1_aux = []
        ninho_2_aux = []
        for j in range(int(len(ninho_x[i]) / 2)):
            ninho_1_aux.append(ninho_x[i][j])
            ninho_2_aux.append(ninho_x[i][j + 10])
        ninho_1.append(ninho_1_aux)
        ninho_2.append(ninho_2_aux)

    for i in range(len(ninho_x)):
        posicoes_1_aux = []
        posicoes_2_aux = []
        for j in range(int(len(ninho_x[0])/2)):
            indice_1 = ninho_1[i].index(max(ninho_1[i]))
            indice_2 = ninho_2[i].index(max(ninho_2[i]))
            posicoes_1_aux.append(indice_1 + 1)
            ninho_1[i][indice_1] = -1e+9
            posicoes_2_aux.append(indice_2 + 11)
            ninho_2[i][indice_2] = -1e+9
        posicoes_1.append(posicoes_1_aux)
        posicoes_2.append(posicoes_2_aux)

    for i in range(len(ninho_x)):
        for j in range(int(len(ninho_x[0])/2)):
            posicoes[i][j] = posicoes_1[i][j]
            posicoes[i][j + 10] = posicoes_2[i][j]

    # print('Entrada_1: {}'.format(posicoes_1))
    # print('Entrada 2: {}'.format(posicoes_2))
    # print('Criando Entradas... entradas = {}\n'.format(posicoes))
    return posicoes
###################################################################################################

###################################################################################################
# Função objetivo
def f_obj(ninho_obj):
    fit = np.ones((n, 1))
    entradas = criando_entrada(ninho_obj)
    for i in range(len(entradas)):
        arq = open('Jessica_Saida.dat', 'w')
        for j in range(len(entradas[i])):
            arq.write('{}\n'.format(str(int(entradas[i][j]))))
        arq.close()
        try:
            os.remove('Jessica_Pico_Boro.dat')
        except:
            pass
        os.startfile('Jessica.exe')
        while True:
            try:
                arq = open('Jessica_Pico_Boro.dat')
                content = arq.read().split("\n")
                arq.close()
                if len(content) >= 2:
                    break
            except:
                time.sleep(1.5)

        pico, boro = lendo_boro_pico()
        if pico <= 1.395:
            fit_aux = 1.0 / boro
        else:
            fit_aux = pico
        fit[i][0] = fit_aux
        # print('Entradas: {}'.format(entradas))
        # print('Fitness: {}'.format(fit))
    return fit, entradas
###################################################################################################

###################################################################################################
# Atualizando os melhores ninhos

melhor_fitness = 10e9

def melhor_ninho(ninho, ninho_novo, fitness):
    fitness_nova, entradas = f_obj(ninho_novo)

    fit_bool = fitness_nova <= fitness # Se a fitness_nova <= fitness anterior, fit_bool = True

    fitness[fit_bool] = fitness_nova[fit_bool] # Se o valor for menor, substitui, se não, mantem valor antigo

    global melhor_fitness


    if fitness.min() < melhor_fitness:
        # print(melhor_fitness)
        melhor_fitness = fitness.min()
        global melhor_posicao
        melhor_posicao = entradas[fitness.argmin()]

    ninho[fit_bool.flatten()] = ninho_novo[fit_bool.flatten()]
    # Transforma de matriz em vetor e susbtitui de acordo com o fit_bool
    return {'min': fitness.min(), 'melhor': fitness.argmin(), 'ninho': ninho, 'fitness': fitness, 'entradas': melhor_posicao}
###################################################################################################

###################################################################################################
def voo_levy(ninho, melhor, li, ls):
    beta = 1.5 # Parâmetro da distribuição de Levy
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2.**((beta-1)/2.)))**(1/beta)
    # Desvio padrão da distribuição do parâmetro u (do parâmetro v é 1)

    # Parâmetros da distribuição de Levy dados através de uma distribuição normal e média zero
    u = np.random.normal(0., sigma, np.shape(ninho))
    v = np.random.normal(0., 1., np.shape(ninho))

    passo = u/((abs(v))**(1/beta)) # Tamanho do salto do voo de Levy
    alfa = 0.75 # alfa = random.randint()
    # Parâmetro da distribuição de Levy 
    ninho = ninho + alfa*passo # Atualizando os ninhos

    # Testando os limites
    bli = ninho < li # Se o ninho for menor que o limite inferior bli = True
    bls = ninho > ls # Se o ninho for maior que o limite superior bls = True

    ninho[bli] = li[bli] # Se a posição correspondente a bli no ninho for True substitui em li
    ninho[bls] = ls[bls] # Se a posição correspondente a bls no ninho for True substitui em ls

    return ninho
###################################################################################################

###################################################################################################
# Voo Aleatório = Abandono dos Ninhos
def voo_aleatorio(ninho, li, ls, pa):
    K = np.random.rand(np.shape(ninho)[0], np.shape(ninho)[1]) < pa # Se a matriz aleatória for menor que pa, K = True
    ninho_novo = ninho + K *(np.random.permutation(ninho) - np.random.permutation(ninho))*np.random.rand(np.shape(ninho)[0], np.shape(ninho)[1])
    # Atualizando ninhos novos depois do abandono

    bli = ninho < li # Se ninho < li então bli = True
    bls = ninho > ls # Se ninho > ls então bls = True

    ninho_novo[bli] = li[bli] # Se a posição correspondente a bli no ninho_novo for True substitui em li
    ninho_novo[bls] = ls[bls] # Se a posição correspondente a bli no ninho_novo for True substitui em lio

    return ninho_novo
###################################################################################################

# Criando o loop
r = melhor_ninho(ninho, ninho, fitness)
N = 0
# print('Dicionário: {}\n\n'.format(r))

while(r['min'] > T):
    lista = []
    for i in range(len(r['fitness'])):
        lista.append(r['fitness'][i][0])
    print('\t\33[35mFitness: {}\33[m'.format(lista))

    r_voo = voo_levy(r['ninho'], r['melhor'], li, ls)

    r = melhor_ninho(r['ninho'], r_voo, r['fitness'])

    r_voo_aleatorio = voo_aleatorio(r['ninho'], li, ls, pa)

    r = melhor_ninho(r['ninho'], r_voo_aleatorio, r['fitness'])
    print('Ciclo {} de {} possíveis\tMenor Valor de Fitness = \33[36m{}\33[m'.format(N, max_geracoes, r['min']))
    N += 1
    if N > max_geracoes:
        break

arq = open('Jessica_Saida.dat', 'w')
for j in range(len(r['entradas'])):
    arq.write('{}\n'.format(str(int(r['entradas'][j]))))
arq.close()
os.startfile('Jessica.exe')
time.sleep(2.5)

t_final = time.ctime()
pico, boro = lendo_boro_pico()
arq_save = open('Jessica_Melhor_configuração_{}_{}.dat'.format(n, max_geracoes), 'w')
arq_save.write('### Programa executado com {} ninhos por geração durante {} gerações\n###\t'
               'Probabilidade de Abandono: {}\n###\tFitness: {}\n###\tInício da Execução: {}\t Final da Execução: {}\n'
               '###\tFator de Pico: {}\tConcentração de Boro: {}\n'.format(n, max_geracoes, pa, r['min'], t_inicio, t_final, pico, boro))
for i in range(len(r['entradas'])):
    arq_save.write('{}\n'.format(int(r['entradas'][i])))
arq_save.close()
print('Parou no loop: {}'.format(N-1))
print('\n\33[35mMelhor resultado: {}\33[m'.format(r['min']))