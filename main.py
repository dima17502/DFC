import re
from copy import copy, deepcopy

'''
   
'''
func_list = []
vector = ""
anf = ""
walsh_coefs = []
hamming_weight = ''
degree = ''
is_balanced = False 
is_affin = False 


def get_func():    
    global vector, anf, walsh_coefs
    print("Выберите способ задания булевой функции: \n1 - вектор значений\n2 - вектор коэффициентов Уолша-Адамара\n3 - В виде алгебраической нормальной формы")
    user_input = input()
    if user_input == "1":
        print("Введите вектор значений булевой функции:")
        vector = input() 
        anf = vector_to_anf(vector)
        walsh_coefs = vector_to_walsh(vector)
    elif user_input == "2":
        print("Введите вектор коэффициентов Уолша-Адамара:")
        walsh_coefs = list(map(int, input().split()))
        vector = walsh_to_vector(walsh_coefs) 
        print(vector)
        anf = vector_to_anf(vector)
        print(anf)
        
    elif user_input == "3":
        print("Введите запись булевой функции (x1,x2,...xn) в виде полинома Жегалкина. Например: 1+x1x2+x1x2x3 - запись для f(x1,x2,x3), '+' соответсвует операции XOR, знак умножения пропускается")
        anf = input()
        vector = anf_to_vector(anf)
        walsh_coefs = vector_to_walsh(vector)
    return vector

def get_adamar_matrix(n):
    '''
        Cложность алгоритма C*2^(2n)
    '''

    if n == 1:
        return [[1, 1],[1, -1]]
    else:
        k = 1
        res = [[1,1],[1,-1]]
        while k < n:
            hn = deepcopy(res)
            for i in range(len(hn)):
                for j in range(len(hn)):
                    res[i].append(hn[i][j])
            copy_res = deepcopy(res)
            for i in range(len(hn)):
                res.append(copy_res[i])
            for i in range(len(res)//2, len(res)):
                for j in range(len(res)//2, len(res)):
                    res[i][j] = -1*res[i][j]
            k += 1
        return res

def count_mutual_correlation(v1, v2, epsilon):
    # Сложность O(n*2^n)
    total = 0
    for i in range(len(v1)): # 2^n

        c = n_to_binary(i)
        c = '0'* (len(epsilon) - len(c)) + c 
        beta_eps = ''
        for j in range(len(c)): # n 
            beta_eps += str(int(c[j])^int(epsilon[j]))
        total += (-1)**(int(v1[i])^int(v2[bin_to_n(beta_eps)]))
    return total

def bin_to_n(vector):
    # Сложность: O(C*2^n)
    total = 0
    p = 0
    for c in vector[::-1]:  # 2^n
        total += int(c)*2**(p)
        p += 1
    return total

def vector_to_walsh(vector):
    '''
        сложность алгоритма ~ 2*2^(2n)
    '''
    h_n = get_adamar_matrix(is_deg_of_two(len(vector)))
    col = []
    for c in vector:
        col.append((-1)**int(c))
    w_coefs = []
    for j in range(len(h_n)):
        total = 0
        for i in range(len(h_n)):
            total += col[i]*h_n[i][j]
        w_coefs.append(total)
    return w_coefs

def walsh_to_vector(walsh_coefs):
    '''
      Сложность: ~ O(2^(2n+2))
    '''

    h_n = get_adamar_matrix(is_deg_of_two(len(walsh_coefs)))        # ~ 2^(2n+1)
    values = []
    for j in range(len(h_n)):           # 2^n
        total = 0
        for i in range(len(h_n)):       # * 2^n
            total += h_n[i][j]* walsh_coefs[i]
        values.append(total)
    vector = []
    for v in values:        # 2^n
        if (v // len(h_n) ) == -1:
            vector.append(1)
        else:
            vector.append(0)

    return vector

def anf_to_vector(anf):
    #  Сложность: O(2^(2n))
    mat = []
    mat.append(anf)
    for i in range(1, len(anf)):        # 2^n
        mat.append([])
        for j in range(len(anf) - i):   # 2^(n-1)
            mat[i].append(int(mat[i-1][j])^int(mat[i-1][j+1]))
    vector = ''
    for v in mat:   # 2^n
        vector += str(v[0])
    return vector

def count_properties(vector):
    properties = {}
    w = 0
    for c in vector:
        w += int(c)
    properties["вес"] = w
    if w == len(vector) // 2:
        properties["is_balanced"] = True 
    else:
        properties["is_balanced"] = False 

    degree = 0 
    anf  = vector_to_anf(vector)
    #print(anf)
    for i in range(len(anf)):
        if anf[i] == "1":
            units = n_to_binary(i).count('1')
            if units > degree:
                degree = units
    properties['degree'] = degree
    properties['is_affin'] = (degree == 1)

    return properties

def hamming_distance(v1, v2):
    # Сложность O(2^(2n))
    count = 0
    times = 1 
    lower_vector = v2
    longer_vector = v1 
    if is_deg_of_two(v1) and is_deg_of_two(v2):
        if (len(v2) > len(v1)):
            times = len(v2) // len(v1)
            lower_vector = v1
            longer_vector = v2
        else: 
            times = len(v1) // len(v2)
    for i in range(len(lower_vector)):              # 2^n
        for c in longer_vector[i*times:(i+1)*times]:    # 2^n
            if lower_vector[i] != c:
                count += 1
    return count

def n_m_mapping():
    # Сложность O(m * 2^n)
    global func_list
    print("Введите количество булевых функций в отображении:")
    m = int(input())
    for i in range(m):
        func_list.append(get_func())
    print(func_list)
    res = []
    for i in range(len(func_list[0])):
        res.append([])
        for j in range(m):
            res[i].append(func_list[j][i])
    return res 


def vector_to_anf(vector):
    #  Cложность: 2^(2n)
    global anf
    t_m = [vector]
    col = 1
    for i in range(len(vector) - 1,0, -1):      # 2^n
        t_m.append([])
        for j in range(i): # 2^(n-1)
            t_m[col].append(str(int(t_m[col-1][j])^int(t_m[col-1][j+1])))
        col += 1
    res = ''
    for c in t_m: # 2^n
        res += str(c[0])
    anf = res 
    s = ''
    for i in range(len(anf)):       # 2^n
        if anf[i] == '1':
            b = n_to_binary(i)
            zeros = is_deg_of_two(len(vector)) - len(b)
            b = '0'*zeros + b
            if b.find('1') == -1:
                s += '1'
            #print(b)
            for ind in range(len(b)):   # n
                if b[ind] == "1":
                    s += ("x" + str(ind+1))
            s += '⊕ '
    s = s[:-2]
    print(s)
    return res 

def is_deg_of_two(n):
    # Сложность: log2_n
    i = 1 
    cur_num = i
    k = 0
    while cur_num < n:
        cur_num *= 2
        k += 1
    if cur_num != n:
        return False 
    else:
        return k
    

def n_to_binary(n):
    # Сложность log2_n
    res = ''
    if n == 0:
        return '0'
    elif n == 1:
        return '1'
    else:
        while n != 1:
            res += str(n%2)
            n //= 2 
        res += '1'
        res = res[::-1]
        return res 

def main():
    v = get_func()
    print(vector_to_anf(v))
    #anf_to_vector(v)
    #print(get_adamar_matrix(4))
    #w = vector_to_walsh(v)
    #print(walsh_to_vector(w))
    #print(count_properties(v))
    #print(n_m_mapping())
    #print(bin_to_n('111'))
    #print(vector_to_anf(v))
    #print(count_mutual_correlation("10101010","10100101","000"))
    #print(n_to_binary(255))


if __name__ == '__main__':
    main()
