import matplotlib.pyplot as plt
import numpy as np
import func as fu
import matplotlib.ticker as ticker

#===============================================================================================================

nums = fu.Reading('data.txt')

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
theta = diff(w, x)

#print(eq4, end='\n\n\n')

#print(w, end='\n\n\n')

#--------------------------------------------------CONDITIONS-----------------------------------------------------

cond1 = w.subs(x = x_r1) == 0
cond2 = w.subs(x = x_r2) == 0
cond3 = w.subs(x = x_r3) == 0
cond4 = w.subs(x = x_z)  == 0
cond5 = w.subs(x = x_k)  == w_0
cond6 = r_k == -k * w_0                                                        # знак и-за того что при отрицательном прогибе сила должна вверх
cond7 = theta.subs(x = x_z) == 0

cond8 = r_k + r1 + r2 + r3 + r_z + (x_q2 - x_q1) * (q1 + q2)/2 == 0
cond9 = m_0 + r1 * x_r1 + r2 * x_r2 + r3 * x_r3 + p * x_p + m + m_z + r_z * x_z + q1 * (x_q2^2 - x_q1^2) / 2 + k_tg * (x_q2^3 - x_q1^3) / 3 == 0   

#--------------------------------------------------SOLUTION--------------------------------------------------------

solution = solve([cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8, cond9], r1, r2, r3, r_z, m_z, m_0, w_0, theta_0, r_k)
# проблема система не дорешивает и выражает все относительно r1

#print(solution[0], end='\n\n\n')

r1 = solve([solution[0][0]], r1)[0].rhs()                                      # решаю отдельно первое уравнение

values = {eq.lhs() : eq.rhs().subs(r1=r1) for eq in solution[0]}

#print(values, end='\n\n\n')

r2, r3, r_z, m_z, m_0, w_0, theta_0, r_k = values[r2], values[r3], values[r_z], values[m_z], values[m_0], values[w_0], values[theta_0], values[r_k]

#--------------------------------------------------ПРОГИБ------------------------------------------------------------

w = w.subs(r1=r1, r2=r2, r3=r3, r_z=r_z, m_z=m_z, m_0=m_0, w_0=w_0, theta_0=theta_0, r_k=r_k)
theta = diff(w,x)
M = diff(theta,x)

print('w = ', w, end='\n\n\n')
print('theta = ', theta, end='\n\n\n')
print('M = ', M, end='\n\n\n')

#--------------------------------------------------GRAPHICS-----------------------------------------------------------

x_mesh = np.linspace(0,25, 100)

w_lambda = lambda x_ : float(w.subs(x = x_))
theta_lambda = lambda x_: float(theta.subs(x = x_))
M_lambda = lambda x_: float(M.subs(x = x_))

w_max = max([w_lambda(x_) for x_ in x_mesh]);        norm_1 = w_max / (abs(q1) + abs(q2) + 1e-7) 

fig, axes = plt.subplots(3, figsize=(10, 8))

for ax in axes:
    ax.plot([0,25], [0,0], color='purple')
    ax.plot([x_q1, x_q2], [q1 * norm_1, q2 * norm_1], color='blue')
    ax.fill_between([x_q1, x_q2], [q1 * norm_1, q2 * norm_1], color='lightblue', alpha=0.5)

    ax.scatter([x_r1, x_r2, x_r3], [0,0,0], s=50, c='red', marker='^', label='опоры')
    ax.scatter([x_m], [0], s=50, c='green', marker='D', label='момент')
    ax.plot([x_z, x_z], [0.3*w_max, -0.3*w_max], c='Brown', label='заделка', lw=4)
    ax.scatter([x_k], [0], s=50, c='Orange', marker='s', label='пружина')

    ax.set_xlim(0 ,25)
    
    for spine in ax.spines.values():
        spine.set_linewidth(0.2)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.grid(which='major', linewidth=1)
    ax.legend()


axes[0].plot(x_mesh, [w_lambda(x_) for x_ in x_mesh]);          axes[0].set_title('Прогиб')
axes[1].plot(x_mesh, [theta_lambda(x_) for x_ in x_mesh]);      axes[1].set_title('Углы')
#axes[2].plot(x_mesh, [M_lambda(x_) for x_ in x_mesh]);          axes[2].set_title('Момент')


plt.tight_layout()
plt.savefig('q.png')