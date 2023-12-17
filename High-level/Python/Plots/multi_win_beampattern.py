import numpy as np
import matplotlib.pyplot as plt
# import scienceplots
from scipy.signal import windows
# plt.style.use(['science','no-latex'])

nr_sensor_elem = 16 #Number of sensor elements
d = 0.5 # sensor element distance (pitch) d = 0.5 => d=lambda/2 defautl
theta_s = [90, 60, 0]
theta_rad = np.linspace(-np.pi, np.pi, 10000)  
def calc_beam_pattern(window_func, theta_rad, d):
    k = 2 * np.pi / 1  
    array_factor = np.zeros_like(theta_rad, dtype=complex)
    for n in range(nr_sensor_elem):
        # array_factor += window_func[n] * np.exp(1j * n * k * d * np.sin(theta_rad))
        array_factor += window_func[n] * np.exp(1j * n * k * d * np.cos(theta_rad))
    return 20 * np.log10(np.abs(array_factor) / np.abs(array_factor).max())

weight_hamming = windows.hamming(nr_sensor_elem)
weight_uniform = np.ones(nr_sensor_elem)
weight_cos = np.cos(np.pi * np.arange(-nr_sensor_elem/2 + 0.5, nr_sensor_elem/2 + 0.5) / nr_sensor_elem) 
weight_blackman = windows.blackman(nr_sensor_elem) 
weight_kaiser = windows.kaiser(nr_sensor_elem, beta=14) 
# weight_hamming = 0.54347 * (0.54 + 0.46 * np.cos(2 * np.pi * np.arange(-nr_sensor_elem/2 + 0.5, nr_sensor_elem/2 + 0.5) / nr_sensor_elem))
weight_barlett = windows.bartlett(nr_sensor_elem)


plt.figure(figsize=(10, 6))
plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_hamming, theta_rad, d), linewidth ='1.5', label='Hamming', linestyle='solid', color='black')
plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_uniform, theta_rad, d), linewidth ='1.5' , label='Uniform', linestyle='solid', color='blue')
# # plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_cos, theta_rad, d), linewidth ='1.5' ,label='Cosine apodization function', linestyle='solid', color='red')
plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_blackman, theta_rad, d), linewidth ='1.5' ,label='Blackman', linestyle='solid', color='red')
plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_kaiser, theta_rad, d), linewidth ='1.5' ,label='Kaiser', linestyle='solid', color='darkgreen')
# plt.plot(theta_rad * 180 / np.pi, calc_beam_pattern(weight_barlett, theta_rad, d), linewidth ='1.5' ,label='Barlett', linestyle='solid', color='green')


plt.xlim(0, 180)
plt.ylim(-80, 0)
plt.xlabel('Angle deg[$\degree$]')
plt.ylabel('Beam Pattern gain [dB]') 
plt.title('Beam Pattern for different apodization/weighting function. N=16 $d = \lambda / 2$')
plt.legend()
plt.grid(True)
plt.savefig("multw_bp.pdf", format="pdf", bbox_inches="tight")
plt.show()
