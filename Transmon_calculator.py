import math as m
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
import numpy as np

e=1.60217662e-19 # elementary charge in C
h=6.62607015e-34 # Planck constant in J s
hbar=1.0545718e-34 # reduced Planck constant in J s
Delta_Al=170e-6*e # superconducting energy gap of Al in J
Phi0=h/(2*e) # magnetic flux quantum in Wb
RRR = 1.25 # Residual-resistance ratio
Z0 = 50
Zcpw = 47.8


def calculate():
    try:
        # Get input values from the GUI
        L_J = float(L_J_entry.get()) * 1e-9
        C_d = float(C_d_entry.get()) * 1e-15
        C_g = float(C_g_entry.get()) * 1e-15
        C_q = float(C_q_entry.get()) * 1e-15
        w_q = float(w_q_entry.get()) * 2 * m.pi * 1e9
        w_r = float(w_r_entry.get()) * 2 * m.pi * 1e9
        chi = float(chi_entry.get()) * 2 * m.pi * 1e6
        alpha = float(alpha_entry.get()) * 2 * m.pi * 1e6
        C_c = float(C_c_entry.get()) * 1e-15
        C_r = float(C_r_entry.get()) * 1e-15
        C_f = float(C_f_entry.get()) * 1e-15
        C_in = float(C_in_entry.get()) * 1e-15
        w_f = float(w_f_entry.get()) * 2 * m.pi * 1e9
        
        # Calculate output variables
        Delta = w_q - w_r
        g = m.sqrt(chi / (1/Delta - 1/(Delta+alpha)))
        I_c = 1/L_J * Phi0 / 2 / m.pi
        R_n = m.pi/2 * Delta_Al / e / I_c
        R_RT = RRR * R_n
        kappa_f = Z0 * w_f**2 * C_in**2 / C_f
        Q_f = w_f / kappa_f
        G = C_c / 2 * m.sqrt(w_r*w_f / (C_r+C_c+C_g) / (C_f+C_c+C_in))
        kappa_q = 4*abs(G)**2 / kappa_f / (1+(2*(w_f-w_q)/kappa_f)**2)
        T_Purcell = 2*m.pi / kappa_q
        kappa_eff = 4*abs(G)**2 / kappa_f / (1+(2*(w_f-w_r)/kappa_f)**2)
        Omega_Rabi = C_d/(C_d+C_g+C_q)*m.sqrt(hbar / 2 / m.sqrt(L_J/(C_d+C_g+C_q)))/hbar

        # Display the results in a message box
        inputs = f"Input variables (from Qiskit simulation)\n"\
                f"L_J [nH]: {round(float(L_J_entry.get()),3)}\n"\
                f"C_d [fF]: {round(float(C_d_entry.get()),3)}\n"\
                f"C_g [fF]: {round(float(C_g_entry.get()),3)}\n\n"\
                f"C_q [fF]: {round(float(C_q_entry.get()),3)}\n"\
                f"w_q [GHz]: {round(float(w_q_entry.get()),3)}\n"\
                f"w_r [GHz]: {round(float(w_r_entry.get()),3)}\n"\
                f"chi [MHz]: {round(float(chi_entry.get()),3)}\n"\
                f"alpha [MHz]: {round(float(alpha_entry.get()),3)}\n\n"\
                f"C_c [fF]: {round(float(C_c_entry.get()),3)}\n"\
                f"C_r [fF]: {round(float(C_r_entry.get()),3)}\n"\
                f"C_f [fF]: {round(float(C_f_entry.get()),3)}\n\n"\
                f"C_in [fF]: {round(float(C_in_entry.get()),3)}\n\n"\
                f"w_f [GHz]: {round(float(w_f_entry.get()),3)}\n"\
                

        result = f"Delta/2pi [GHz]: {round(Delta/(2*m.pi)/1e9,6)}\n" \
                 f"g/2pi [MHz]: {round(g/(2*m.pi*1e6),6)}\n" \
                 f"G/2pi [MHz]: {round(G/(2*m.pi)/1e6,6)}\n\n" \
                 f"Q_f: {round(Q_f,6)}\n\n" \
                 f"kappa_f/2pi [MHz]: {round(kappa_f/(2*m.pi)/1e6,6)}\n" \
                 f"kappa_q/2pi [MHz]: {round(kappa_q/(2*m.pi)/1e6,6)}\n" \
                 f"kappa_eff/2pi [MHz]: {round(kappa_eff/(2*m.pi)/1e6,6)}\n\n" \
                 f"T_Purcell [us]: {round(T_Purcell*1e6,6)}\n\n" \
                 f"I_c[nA]: {round(I_c*1e9,6)}\n" \
                 f"R_RT[kOhm]: {round(R_RT/1e3,6)}\n" \
                 f"R_n [kOhm]: {round(R_n/1e3,6)}\n\n" \
                 f"Omega_Rabi [MHz] \n for -130[dBm]~100[nV] drive:\n {round(Omega_Rabi/1e6*1e-7,6)}\n\n" \
                 f"Discremination chi / kappa_eff: {round(chi/kappa_eff,3)}\n\n"
        
        fig, ax = plt.subplots(figsize=(10,5))   
        
        ax.set_xlabel('Delta/kappa_eff')
        ax.set_ylabel('g/kappa_eff')
        ax.set_title("Qubit regime")

        # Calculate the data point
        Delta_kappa_eff = abs(Delta / kappa_eff)
        g_kappa_eff = abs(g / kappa_eff)

        # Set the x and y ranges
        ax.set_xlim(0, 1500)
        ax.set_ylim(0, 150)

        # Plot the data point
        ax.scatter(Delta_kappa_eff, g_kappa_eff)

        ax.text(x=1550, y=20, s=inputs)
        ax.text(x=2250, y=0, s=result)
        
        ax.text(x=1550, y=10, s='Delta/kappa_eff :  '+str(Delta_kappa_eff))
        ax.text(x=1550, y=0, s='g/kappa_eff :  '+str(g_kappa_eff))
        

        x_pts = np.linspace(0,1500,151)
        guide1 = np.sqrt(x_pts)
        guide2 = x_pts / 10
        guide3 = x_pts
        guide4 = x_pts**(3/4)
        
        ax.plot(x_pts,guide1,'r')
        ax.plot(x_pts,guide2,'g')
        ax.plot(x_pts,guide3,'b')
        ax.plot(x_pts,guide4,'k')

        plt.show()

        # messagebox.showinfo("Results", result)


    except ValueError:
        messagebox.showerror("Error", "Invalid input. Please enter numeric values.")

# Create the GUI window
window = tk.Tk()
window.title("Transmon Calculator")

# Create input labels and entry fields
L_J_label = tk.Label(window, text="L_J [nH]:")
L_J_label.grid(row=1, column=0)
L_J_entry = tk.Entry(window)
L_J_entry.grid(row=1, column=1)






C_d_label = tk.Label(window, text="C_d [fF]:")
C_d_label.grid(row=3, column=0)
C_d_entry = tk.Entry(window)
C_d_entry.grid(row=3, column=1)

C_g_label = tk.Label(window, text="C_g [fF]:")
C_g_label.grid(row=4, column=0)
C_g_entry = tk.Entry(window)
C_g_entry.grid(row=4, column=1)

C_q_label = tk.Label(window, text="C_q [fF]:")
C_q_label.grid(row=5, column=0)
C_q_entry = tk.Entry(window)
C_q_entry.grid(row=5, column=1)

w_q_label = tk.Label(window, text="w_q [GHz]:")
w_q_label.grid(row=7, column=0)
w_q_entry = tk.Entry(window)
w_q_entry.grid(row=7, column=1)

w_r_label = tk.Label(window, text="w_r [GHz]:")
w_r_label.grid(row=8, column=0)
w_r_entry = tk.Entry(window)
w_r_entry.grid(row=8, column=1)

chi_label = tk.Label(window, text="chi [MHz]:")
chi_label.grid(row=9, column=0)
chi_entry = tk.Entry(window)
chi_entry.grid(row=9, column=1)

alpha_label = tk.Label(window, text="alpha [MHz]:")
alpha_label.grid(row=10, column=0)
alpha_entry = tk.Entry(window)
alpha_entry.grid(row=10, column=1)

C_c_label = tk.Label(window, text="C_c [fF]:")
C_c_label.grid(row=3, column=3)
C_c_entry = tk.Entry(window)
C_c_entry.grid(row=3, column=4)

C_r_label = tk.Label(window, text="C_r [fF]:")
C_r_label.grid(row=4, column=3)
C_r_entry = tk.Entry(window)
C_r_entry.grid(row=4, column=4)

C_f_label = tk.Label(window, text="C_f [fF]:")
C_f_label.grid(row=5, column=3)
C_f_entry = tk.Entry(window)
C_f_entry.grid(row=5, column=4)

C_in_label = tk.Label(window, text="C_in [fF]:")
C_in_label.grid(row=6, column=3)
C_in_entry = tk.Entry(window)
C_in_entry.grid(row=6, column=4)

w_f_label = tk.Label(window, text="w_f [GHz]:")
w_f_label.grid(row=8, column=3)
w_f_entry = tk.Entry(window)
w_f_entry.grid(row=8, column=4)

# Create the calculate button
calculate_button = tk.Button(window, text="Calculate", command=calculate)
calculate_button.grid(row=15, column=0, columnspan=2)

# Set preset values for inputs

A1 = tk.Label(window, text="*----Junction setup----*")
A1.grid(row=0, column=1)
L_J_entry.insert(0, "9")

A2 = tk.Label(window, text="*----C-sim Q----*")
A2.grid(row=2, column=1)
C_d_entry.insert(0, "2.01854")
C_g_entry.insert(0, "22.6872")
C_q_entry.insert(0, "122.913")

A3 = tk.Label(window, text="*----pyEPR----*")
A3.grid(row=6, column=1)
w_q_entry.insert(0, "4.66608")
w_r_entry.insert(0, "5.55215")
chi_entry.insert(0, "-5.31")
alpha_entry.insert(0, "-171.77")

A4 = tk.Label(window, text="*----C-sim R-PF-TL----*")
A4.grid(row=2, column=4)
C_c_entry.insert(0, "4.39624")
C_r_entry.insert(0, "1873.67")
C_f_entry.insert(0, "1874.82")
C_in_entry.insert(0, "79.0495")

A5 = tk.Label(window, text="*----S21 analysis----*")
A5.grid(row=7, column=4)
w_f_entry.insert(0, "5.621")

# Start the GUI event loop
window.mainloop()

