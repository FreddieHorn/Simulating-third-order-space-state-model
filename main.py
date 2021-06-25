
import matplotlib
import numpy as np
import tkinter as tk
import tkinter.messagebox
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import tkinter.font as font
import math
import cmath

#funkcja obliczajaca faze i amo
def solving_mag_phase(a0, a1, a2, b0, b1, b2, b3, h):
    G_phase = []
    G_mag = []
    omega = []
    w = np.arange(0, 20 + h, h)
    for i in w:
        if i == 0 and b0 == 0 or (a0-a2*math.pow(i, 2) == 0 and a1 * i - math.pow(i, 3) == 0):
            continue
        else:
            Gjw = (complex(b0-b2*math.pow(i, 2), b1*i - b3 * math.pow(i, 3))/(complex(a0-a2*math.pow(i, 2), a1 * i - math.pow(i, 3))))
            G_phase.append(cmath.phase(Gjw)*(180/math.pi))
            G_mag.append(20*cmath.log10(abs(Gjw)))
            omega.append(i)

    return G_phase, G_mag, omega

#funkcja wykonująca główny algorytm
def solving(A, B, C, D, u, t, h):
    x = np.zeros((3, 1))
    xdt = []
    y = np.zeros(len(t))

    for i in range(len(t)):
        Ax = np.dot(A, x)
        Bu = np.dot(B, u[i])
        Cx = np.dot(C, x)
        Du = np.dot(D, u[i])
        xdt = np.add(Ax, Bu)
        xdt = xdt*h
        xdt = np.add(x, xdt)
        x = xdt
        y[i] = np.add(Cx, Du)

    return y

#w tej funkcji przypisujemy liczby wpisane w interfejsie do zmiennych oraz pokazujemy wykresy
def show_plots():

    if T_t.get() <= 0:
        tk.messagebox.showwarning(title="Blad", message="Czas trwania symulacji musi byc wiekszy od zera")
    elif h_t.get() <=0:
        tk.messagebox.showwarning(title="Blad", message="Krok symulacji musi być większy od zera")
    else:
        if a0_t.get() < 0 or a1_t.get() < 0 or a2_t.get() < 0 or a0_t.get()-a1_t.get()*a2_t.get()>=0:
            tk.messagebox.showwarning(title="Blad",
                                      message="Układ niestabilny!")
        a0 = a0_t.get()
        a1 = a1_t.get()
        a2 = a2_t.get()
        b0 = b0_t.get()
        b1 = b1_t.get()
        b2 = b2_t.get()
        b3 = b3_t.get()

        T = T_t.get()  # czas trwania symulacji
        h = h_t.get()  # krok symulaji
        M = M_t.get()  # amplituda sygnału wejściowego
        L = L_t.get()

        w = L/T
        t = np.arange(0, T + h, h)

        A = np.matrix([[0, 0, -a0],
                       [1, 0, -a1],
                       [0, 1, -a2]])

        B = np.matrix([[b0 - a0 * b3],
                       [b1 - a1 * b3],
                       [b2 - a2 * b3]])

        C = np.matrix([0, 0, 1])

        D = np.matrix([b3])

        u = M * np.sign(np.sin(t * 2 * np.pi*w))
        u1 = M * np.ones(len(t))
        u2 = M * np.sin(t * 2 * np.pi*w)

        if v.get() == "1":
            wejscie = u
        elif v.get() == "2":
            wejscie = u1
        else:
            wejscie = u2

        matplotlib.pyplot.close()

        plt.figure()
        plt.subplot(221)
        plt.plot(t, solving(A, B, C, D, wejscie, t, h))
        plt.title('Wyjście')
        plt.xlabel('Czas')
        plt.ylabel('Amplituda')
        plt.grid(True)

        plt.subplot(222)
        plt.plot(t, wejscie)
        plt.title('Wejście sygnału')
        plt.xlabel('Czas')
        plt.ylabel('Amplituda')
        plt.grid(True)

        plt.subplot(223)
        plt.plot(solving_mag_phase(a0,a1,a2,b0,b1,b2,b3, h)[2], solving_mag_phase(a0,a1,a2,b0,b1,b2,b3,h)[0])
        plt.xscale('log')
        plt.title('Wykres fazowy')
        plt.xlabel('Omega')
        plt.ylabel('Faza [°]')
        plt.grid(True)

        plt.subplot(224)
        plt.plot(solving_mag_phase(a0,a1,a2,b0,b1,b2,b3,h)[2], solving_mag_phase(a0,a1,a2,b0,b1,b2,b3,h)[1])
        plt.xscale('log')
        plt.title('Wykres amplitudowy')
        plt.xlabel('Omega')
        plt.ylabel('Amplituda [dB]')
        plt.grid(True)
        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.35,
                            wspace=0.35)
        plt.show()


#interfejs
window = tk.Tk()
window.geometry("1000x550")


tk.Label(window, text="a0:").grid(row=0)
tk.Label(window, text="a1:").grid(row=1)
tk.Label(window, text="a2:").grid(row=2)
tk.Label(window, text="b0:").grid(row=3)
tk.Label(window, text="b1:").grid(row=4)
tk.Label(window, text="b2:").grid(row=5)
tk.Label(window, text="b3:").grid(row=6)


a0_t = tk.DoubleVar()
a1_t = tk.DoubleVar()
a2_t = tk.DoubleVar()
b0_t = tk.DoubleVar()
b1_t = tk.DoubleVar()
b2_t = tk.DoubleVar()
b3_t = tk.DoubleVar()


a0_e = tk.Entry(window, textvariable=a0_t).grid(row=0, column=1)
a1_e = tk.Entry(window, textvariable=a1_t).grid(row=1, column=1)
a2_e = tk.Entry(window, textvariable=a2_t).grid(row=2, column=1)
b0_e = tk.Entry(window, textvariable=b0_t).grid(row=3, column=1)
b1_e = tk.Entry(window, textvariable=b1_t).grid(row=4, column=1)
b2_e = tk.Entry(window, textvariable=b2_t).grid(row=5, column=1)
b3_e = tk.Entry(window, textvariable=b3_t).grid(row=6, column=1)

myFont = font.Font(size=20, weight='bold')

Button1 = tk.Button(window, text="Symuluj", command=show_plots, bg='#0052cc', fg='#ffffff')
Button1['font'] = myFont
Button1.config(width=20, height=5)
Button1.place(x=340, y=355)


Opis_wejscia = tk.Label(window, text="sygnał wejściowy:")
Opis_wejscia.grid(row=0, column=2)

v = tk.StringVar(window, "1")


values = {"Fala prostokątna": "1",
          "Skok jednostkowy": "2",
          "Sygnał sinusoidalny": "3",
          }

for (text, value) in values.items():
    tk.Radiobutton(window, text=text, variable=v,
        value=value).grid(row=int(value), column = 2)

T_t = tk.DoubleVar()
h_t = tk.DoubleVar()
M_t = tk.DoubleVar()
L_t = tk.DoubleVar()

tk.Label(window, text="Czas trwania symulacji: ").grid(row=0, column=3)
tk.Label(window, text="Krok symulacji: ").grid(row=1, column=3)
tk.Label(window, text="Amplituda sygnału wejściowego: ").grid(row=2, column=3)
tk.Label(window, text="Liczba okresów sygnału sinusoidalnego w czasie trwania symulacji: ").grid(row=3, column=3)

T_e = tk.Entry(window, textvariable=T_t)
T_e.grid(row=0, column=4)
h_e = tk.Entry(window, textvariable=h_t)
h_e.grid(row=1, column=4)
m_e = tk.Entry(window, textvariable=M_t)
m_e.grid(row=2, column=4)
L_e = tk.Entry(window, textvariable=L_t)
L_e.grid(row=3, column=4)

#wczytanie zdjęcia
canvas = tk.Canvas(window, width=900, height=200)
canvas.place(x=50, y=160)
img = ImageTk.PhotoImage(Image.open("model.jpg"))
canvas.create_image(0, 10, anchor=tk.NW, image=img)

window.mainloop()
