import sys, random, time, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '..')
import main

N_list = np.array([15, 20, 30, 40, 50, 60, 70, 80, 100])

T_C_list = np.array([2.55, 2.475, 2.41, 2.38, 2.37, 2.35, 2.34, 2.33, 2.29])

T_C_error_list = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.03, 0.03, 0.03, 0.04])

vs = np.linspace(0.5, 4, num = 1000)
R2s, R2s_plus, R2s_minus, cofs, cofs_plus, cofs_minus = [], [], [], [], [], []

def get_R2(x, y, degree):
	

	coeffs = np.polyfit(x, y, degree)

	# r-squared
	p = np.poly1d(coeffs)
	
	# fit values, and mean
	yhat = p(x)                         # or [p(z) for z in x]
	ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
	ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
	R2 = ssreg / sstot

	return R2


for v in vs : 

	#Calculate gradient
	coefs = np.polyfit(N_list**(-1./v), T_C_list, 1)

	coefs_plus = np.polyfit(N_list**(-1./v), T_C_list + T_C_error_list, 1)
	
	coefs_minus = np.polyfit(N_list**(-1./v), T_C_list - T_C_error_list, 1)

	cofs.append(coefs)
	cofs_plus.append(coefs_plus)
	cofs_minus.append(coefs_minus)
	
	#Get regression error from polyfit
	R2 = get_R2(N_list**(-1./v), T_C_list, 1)
	R2_plus = get_R2(N_list**(-1./v), T_C_list + T_C_error_list, 1)
	R2_minus = get_R2(N_list**(-1./v), T_C_list - T_C_error_list, 1)
	
	R2s.append(R2)
	R2s_plus.append(R2_plus)
	R2s_minus.append(R2_minus)

R2 = max(R2s)
R2_max = max(R2s_plus)
R2_min = max(R2s_minus)

v = vs[R2s.index(max(R2s))]
v_max = vs[R2s_plus.index(R2_max)]
v_min = vs[R2s_minus.index(R2_min)]

text_string = "v =  "+ ("%.2f" % v) + " + " 
text_string += ("%.2f" % (v_max-v)) 
text_string +=  " - " + ("%.2f" % (v- v_min))

a, c = np.polyfit(N_list**(-1./v), T_C_list, 1)
a_max, c_max = np.polyfit(N_list**(-1./v_max), T_C_list, 1)
a_min, c_min = np.polyfit(N_list**(-1./v_min), T_C_list, 1)
print()
print("a, a_max, a_min", a, a_max, a_min)
print("c, c_max, c_min", c, c_max, c_min)

print()
print("a =", a, "+", a_min - a, "-", a - a_max)
print("(errors from min v and max v respectively)")
print()
print("Tc(inf) =", c, "+", c_min - c, "-", c- c_max)
print("(errors from min v and max v respectively)")
print()

plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

plt.figure()
plt.plot(vs, R2s)
plt.plot(vs, R2s_plus, '--')
plt.plot(vs, R2s_minus, '--')
plt.legend(["fit", "fit_for_max", "fit_for_min"], prop={ 'size' : 20})

plt.xlabel("Polynomial v", size =22)
plt.ylabel("R-squared fit", size =22)
plt.title("R-square fits for various values of v in N^(-1/v) finite scaling", size =22)

plt.text(vs[int(len(vs)/8)], R2s[int(len(R2s)*0.001)], text_string, size =22 )
plt.rcParams.update({'font.size': 22})

plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

plt.show()

"""
results['polynomial'] = coeffs.tolist()

# r-squared
p = numpy.poly1d(coeffs)
# fit values, and mean
yhat = p(x)                         # or [p(z) for z in x]
ybar = numpy.sum(y)/len(y)          # or sum(y)/len(y)
ssreg = numpy.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
sstot = numpy.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
results['determination'] = ssreg / sstot

return results

"""