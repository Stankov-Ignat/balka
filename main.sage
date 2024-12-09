import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

def reading(s):
    with open(s, 'r') as file:
        nums = []
        for line in file:
            num = line.split('<==')[0].strip()
            try:
                num = float(num)
                nums.append(num)
            except ValueError:
                print('ERROR ЧТЕНИЯ ФАЙЛА')
    return nums

def read_and_save_data(filename, headers):
    data = {header: [] for header in headers}  # Инициализация словаря с указанными ключами

    with open(filename, 'r') as file:
        for line in file:
            values = list(map(float, line.split()))
            for i, value in enumerate(values):
                data[headers[i]].append(value)
    
    return data

def read_numbers_to_list(filename):
    numbers = []  # Инициализируем пустой список для чисел

    with open(filename, 'r') as file:
        for line in file:
            try:
                # Преобразуем строку в число с плавающей точкой и добавляем в список
                numbers.append(float(line.strip()))
            except ValueError:
                # Если строка не является числом, пропускаем её
                continue
    
    return numbers

#===============================================================================================================

nums = reading('data.txt')

x_k = nums[0];     x_q1 = nums[1];    x_q2 = nums[2]
x_r1 = nums[3];    x_r2 = nums[4];    x_r3 = nums[5]
x_m = nums[6];     x_z = nums[7];     x_p = nums[8]

q1 = nums[9];      q2 = nums[10];     p = nums[11];
k = nums[12];      m = nums[13]
E = nums[14];      J = nums[15]
k_tg = (q2 - q1) / (x_q2  - x_q1)

#================================================================================================================

r1, r2, r3, r_z, m_z, m_0, w_0, theta_0, r_k, x = var('r1 r2 r3 r_z m_z m_0 w_0 theta_0 r_k x')

eq1 = m_0 * x^2/2 + r_k * x^3/6
eq2 = heaviside(x - x_q1) * (q1 * (x - x_q1)^4/24 + k_tg * (x - x_q1)^5/120)
eq3 = heaviside(x - x_r1) * (r1 * (x - x_r1)^3/6)
eq4 = heaviside(x - x_p) *  (p * (x - x_p)^3/6)
eq5 = heaviside(x - x_r2) * r2 * (x - x_r2)^3/6
eq6 = heaviside(x - x_q2) * (-q2 * (x - x_q2)^4/24 - k_tg * (x - x_q2)^5/120)
eq7 = heaviside(x - x_r3) * r3 * (x - x_r3)^3/6
eq8 = heaviside(x - x_m) * (m * (x - x_m)^2/2)

w = w_0 + theta_0 * x + 1/(E*J) * (eq1 + eq2 + eq3 + eq4 + eq5 +  eq6 + eq7 + eq8)

theta = diff (w, x)
theta = theta.substitute_function (dirac_delta, lambda x: 0)

M = diff (theta,x)
M = M.substitute_function (dirac_delta, lambda x: 0)

Q = diff (M,x)
Q = Q.substitute_function (dirac_delta, lambda x: 0)

#--------------------------------------------------CONDITIONS-----------------------------------------------------

cond1 = w.subs (x = x_r1) == 0
cond2 = w.subs (x = x_r2) == 0
cond3 = w.subs (x = x_r3) == 0
cond4 = w.subs (x = x_z)  == 0
cond5 = M.subs (x = x_k)  == 0
cond6 = r_k == -k * w_0                                                        # знак и-за того что при отрицательном прогибе сила должна вверх
cond7 = theta.subs (x = x_z) == 0

cond8 = p + r_k + r1 + r2 + r3 + r_z + (x_q2 - x_q1) * (q1 + q2)/2 == 0        # условия равновесия
cond9 = m_0 + r1 * x_r1 + r2 * x_r2 + r3 * x_r3 + p * x_p + m + m_z + r_z * x_z + q1 * (x_q2^2 - x_q1^2) / 2 + k_tg * (x_q2^3 - x_q1^3) / 3 == 0   

#--------------------------------------------------SOLUTION--------------------------------------------------------

solution = solve([cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8, cond9], r1, r2, r3, r_z, m_z, m_0, w_0, theta_0, r_k)

values = {eq.lhs() : eq.rhs() for eq in solution[0]}

for i, (key, value) in enumerate(values.items()):
    if i >= 5:
        break
    print(f"{key}: {value.n():.6f} ", end=" ")
print()

r1, r2, r3, r_z, m_z, m_0, w_0, theta_0, r_k = values[r1],values[r2],values[r3],values[r_z],values[m_z],values[m_0],values[w_0],values[theta_0],values[r_k]

#--------------------------------------------------ПОДСТАНОВКА------------------------------------------------------------

w = w.subs (r1=r1, r2=r2, r3=r3, r_z=r_z, m_z=m_z, m_0=m_0, w_0=w_0, theta_0=theta_0, r_k=r_k)

theta = theta.subs (r1=r1, r2=r2, r3=r3, r_z=r_z, m_z=m_z, m_0=m_0, w_0=w_0, theta_0=theta_0, r_k=r_k)

M = M.subs (r1=r1, r2=r2, r3=r3, r_z=r_z, m_z=m_z, m_0=m_0, w_0=w_0, theta_0=theta_0, r_k=r_k)

Q = Q.subs (r1=r1, r2=r2, r3=r3, r_z=r_z, m_z=m_z, m_0=m_0, w_0=w_0, theta_0=theta_0, r_k=r_k)

#--------------------------------------------------READING CPP RESULT-------------------------------------------------------

headers = ['mke_mesh', 'mke_w', 'mke_theta', 'mke_moment', 'mke_q']

mke_res = read_and_save_data('output.txt', headers)

mesh_orig = read_numbers_to_list ("mesh.txt")

solve = read_numbers_to_list ("solve.txt")

solve_w = solve[0::2];      solve_theta = solve[1::2]

#--------------------------------------------------GRAPHICS-----------------------------------------------------------------

#print ("\nSAGE__GRAFICS")

x_mesh = np.linspace(0,25, 1000)

w_lambda =     lambda x_ : float(w.subs(x = x_))
theta_lambda = lambda x_: float(theta.subs(x = x_))
M_lambda =     lambda x_: float(M.subs(x = x_))
Q_lambda =     lambda x_: float(Q.subs(x = x_))

w_max = max([w_lambda(x_) for x_ in x_mesh]);

norm_q = w_max / (abs(q1) + abs(q2) + 1e-7 + 1e2)                                            # нормализация для распределенной нагрузки

#---------------------------------------------------------------------------------------

fig, axes = plt.subplots (4, 2, figsize = (18, 11))

i = 0
for ax in axes.flatten():

    ax.plot([0,25], [0,0], color='purple')                                             # ось x

    ax.plot([x_q1, x_q2], [0, 0], c='blue', label='q_распр', lw=2)                       # распределенная нагрузка просто на оси

    #ax.plot([x_q1, x_q2], [q1 * norm_q, q2 * norm_q], color='blue')                    # распределенная нагрузка ловкая
    #ax.fill_between([x_q1, x_q2], [q1 * norm_q, q2 * norm_q], color='lightblue', alpha=0.5)

    ax.scatter([x_r1, x_r2, x_r3], [0,0,0], s=50, c='red', marker='^', label='опоры')
    ax.scatter([x_m], [0], s=50, c='green', marker='D', label='момент')
    ax.plot([x_z, x_z], [0.1*w_max, -0.1*w_max], c='Brown', label='заделка', lw=4)
    ax.scatter([x_k], [0], s=50, c='Orange', marker='s', label='пружина')
    ax.scatter([x_p], [0] ,s=50, c='black', marker='^', label='сила')


    ax.set_xlim(0 ,25)
    
    for spine in ax.spines.values():
        spine.set_linewidth(0.2)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.grid(which='major', linewidth=1)
    #ax.legend()

    i += 1

axes[0, 0].legend()
axes[0, 0].plot(x_mesh, [Q_lambda(x_) for x_ in x_mesh]);                              axes[0, 0].set_title('МЕТОД НАЧАЛЬНЫХ ПАРАМЕТРОВ\nСила')
axes[1, 0].plot(x_mesh, [E*J * M_lambda(x_) for x_ in x_mesh]);                        axes[1, 0].set_title('Момент')
axes[2, 0].plot(x_mesh, [theta_lambda(x_) for x_ in x_mesh]);                          axes[2, 0].set_title('Углы')
axes[3, 0].plot(x_mesh, [w_lambda(x_) for x_ in x_mesh]);                              axes[3, 0].set_title('Прогиб')

axes[0, 1].legend()
axes[0, 1].plot (mke_res['mke_mesh'], mke_res['mke_q'], c = 'purple');                 axes[0, 1].set_title('МЕТОД КОНЕЧНЫХ ЭЛЕМЕНТОВ\nCила')
axes[1, 1].plot (mke_res['mke_mesh'], mke_res['mke_moment'], c = 'purple');            axes[1, 1].set_title('Момент')
axes[2, 1].plot (mke_res['mke_mesh'], mke_res['mke_theta'], c = 'purple');             axes[2, 1].set_title('Углы')
axes[3, 1].plot (mke_res['mke_mesh'], mke_res['mke_w'], c = 'purple');                 axes[3, 1].set_title('Прогиб')

# решение из метода бубнова-галеркина
axes[3, 1].scatter (mesh_orig, solve_w, c = 'black', s=20);                 
axes[2, 1].scatter (mesh_orig, solve_theta, c = 'black', s=20);                 

# мкэ на левой стороне
'''
axes[0, 0].plot (mke_res['mke_mesh'], mke_res['mke_q'], c = 'purple');                 
axes[1, 0].plot (mke_res['mke_mesh'], mke_res['mke_moment'], c = 'purple');           
axes[2, 0].plot (mke_res['mke_mesh'], mke_res['mke_theta'], c = 'purple');             
axes[3, 0].plot (mke_res['mke_mesh'], mke_res['mke_w'], c = 'purple');   
'''

plt.tight_layout()
plt.savefig('result.png')

'''
print ("max_w = ", max (mke_res['mke_w']), end = "     ")
print ("max_theta = ", max (mke_res['mke_theta']), end = "     ")
print ("max_moment = ", max (mke_res['mke_moment']), end = "     ")
print ("max_q = ", max (mke_res['mke_q']))
print ("w(0) = ", mke_res['mke_w'][0])
'''