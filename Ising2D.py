# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 12:53:31 2015

@author: Arveiler
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
from time import clock
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk


def initialise():
    """
    Generates a lattice of spins for a temperature of 5.
    """
    L = np.zeros([N,N])+1 # A lattice full of spin-up particules
    for i in range (50):
        L = run_metro(L,5) 
    return L # Lattice after a long relaxation time

def metropolis(L, T):
    """
    Operates the Metropolis algorithm on a lattice.
    L : Lattice.
    T : Temperature.
    """
    i, j = np.random.randint(0,N,size = 2) # Random coordinates in lattice
    inv_spin = L[i][j]  * (-1) # Trial of a new lattice (1 spin reversed)
    accept = False # Boolean to accept or not the new lattice
    n1 = neighbor(L, i, j-1) # 
    n2 = neighbor(L, i-1, j) # Finds the neighbors' spins 
    n3 = neighbor(L, i, j+1) #
    n4 = neighbor(L, i+1, j) #
    var_U = 2 * J * L[i][j] * (n1 + n2 + n3 + n4) # Energy variation due to spin inversion
    if var_U <= 0: # The new lattice is more balanced : it is accepted
        accept = True
    else: # The new lattice is accepted with a Boltzmann factor probability
        p = np.random.random()
        if np.exp(-var_U/(k_b*T)) > p:
            accept = True
    if accept: # The lattice is effectively modified
        L[i][j] = inv_spin
    return L

def neighbor(L, i, j):
    """
    Finds the spin at coordinates (i,j).
    Takes into account the periodicity of the boundary conditions.
    """
    if i < 0:
        i = N-1
    elif i > N-1:
        i = 0
    if j < 0:
        j = N-1
    elif j > N-1:
        j = 0
    return L[i][j]

def run_metro(L,T):
    """
    Runs the Metropolis algorithm (N**2) times so that, statistically, 
    all spins have been picked and reversed once.
    For big-sized lattices, the algorithm is ran fewer times in order
    to save computing time.
    """
    if N<64:
        runs = N*N
    else:
        runs = (N*N)/10
    for i in range (runs):
        L = metropolis(L,T)
    return L

def magnetisation(L):
    """
    Calculates the mean magnetisation by spin of the lattice.
    """
    return abs(float(np.sum(L))/(N**2))

def M_Onsager(var_T):
    """
    Plots the theoretical mean magnetisation of the lattice as function of temperature.
    """
    M_Onsager = []
    for T in var_T:
        if T < Tc:
            M_Onsager.append((1-np.sinh(2.*J/T)**(-4))**(1./8))
        else:
            M_Onsager.append(0)
    plt.plot(var_T,M_Onsager,'k',label = "Onsager's theoretical magnetisation")

def energy(L):
    """
    Calculates the mean energy by spin of the lattice.
    """
    sum_energy = 0
    for i in range(N):
        for j in range(N):
            n1 = neighbor(L, i, j-1)
            n2 = neighbor(L, i-1, j)
            n3 = neighbor(L, i, j+1)
            n4 = neighbor(L, i+1, j)
            sum_energy += L[i][j] * (n1 + n2 + n3 + n4)
    return (-J) * sum_energy / float(2*N**2)

def graphs():
    """
    Displays the evolution of magnetisation and energy of the lattice as functions 
    of temperature as well as Onsager's theoretical magnetisation.
    """
    start = clock()
    H_T = np.linspace(5,2.6,150) # Array of high temperatures
    Tc_T = np.linspace(2.6,1.6,500) # Temperatures around Tc (more points, carrower interval)
    L_T = np.linspace(1.6,0.1,150) # Low temperatures
    var_T = np.concatenate((H_T,Tc_T,L_T)) # Array of decreasing temperatures from 5 to 0
    L = initialise()
    
    # Arrays to fill with the values of magnetisation and energy : 
    M = np.zeros(shape = np.size(var_T))
    U = np.zeros(shape = np.size(var_T))
    i = 0 # Will enable to insert magnetisation and energy values in M and U
    
    # Filling the arrays with values calculated after a long enough relaxation time
    print '\n\n\n\nCalculating...'
    for T in var_T:
        L = run_metro(L,T)
        M[i] = magnetisation(L)
        U[i] = energy(L)
        percentage = i/float(np.size(var_T))*100 # User sees that computing is in progress
        if percentage == 0 or percentage%25 == 0:
            print (percentage),'% done'
        i += 1
    print 'Computing time : ',(clock()-start),'s\n\n\n'
    
    # Creation of the figure and plotting separatly magnetisation and energy
    plt.figure('Magnetisation and Energy',figsize=(15,10))
    
    plt.subplot(2,1,1)
    plt.grid()
    plt.plot(var_T,M,'b',label = 'Magnetisation')
    M_Onsager(var_T) # Plots Onsager's magnetisation
    plt.xlabel('Temperature')
    plt.ylabel('Magnetisation')
    plt.ylim(-.5,1.5)
    # Legends are placed outside the plotting area, on top of it, and are made draggable :
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower center',
           ncol=2, mode="expand", borderaxespad=0.).draggable() 
    
    plt.subplot(2,1,2)
    plt.grid()
    plt.plot(var_T,U,'r',label = 'Energy')
    plt.xlabel('Temperature')
    plt.ylabel('Energy')
    plt.ylim((-2.5,0))
    # Legends are placed outside the plotting area, on top of it, and are made draggable :
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
           ncol=1, mode="expand", borderaxespad=0.).draggable() 
           
    quit_window() # Quits the Tkinter GUI to display the figure
    plt.show()
    
def animated_lattice():
    """
    Displays an animation of the lattice being modified with temperature.
    """
    def update(L):
        """
        Updates the lattice by running the Metropolis algorithm.
        """
        L = run_metro(L,slider.val)
        lattice_image.set_data(L)
        return lattice_image
    
    def generate_lattice():
        yield L # 'yield' instead of 'return' since an iteralble is required in FuncAnimation
    quit_window() # Quits the previous window
    
    # Creation of a figure :
    fig = plt.figure('Ising 2D Animation',figsize=(12,10))
   
   # Creation of a slider in the figure to select the temperature :
    from matplotlib.widgets import Slider
    slax = fig.add_axes((.15,.01,.7,.05))
    slider = Slider(slax,label = 'Temperature',valmin = 0, valmax=5, valinit=5, color = 'gray')
    slax.text(.9,1.1,"Curie's temperature : Tc=2.27")

    # Creation of a plotting area in the figure :
    fig.add_subplot()
    ax = fig.add_axes((.1,.1,.85,.85))

    # Hiding axis text : 
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    # Creation an image of the array-type lattice L :
    L = initialise()
    lattice_image = ax.matshow(L,cmap=cm.gray) 

    # Definition of the animation :
    ani = animation.FuncAnimation(fig, update, generate_lattice, interval = 1)

    plt.show()
    quit_window() # Quits the GUI to display the animation
    return ani
    
    
####################################################################
######################### MAIN #####################################
####################################################################
    
    
# Definition of global variables
global J ; J = 1
global k_b ; k_b = 1
global Tc ; Tc = (2.*J)/((np.log(1+np.sqrt(2)))*k_b)
global N ; N = 0 # The default value is set to zero for later input verification
                 # but the user will choose the real value

set_N = Tk.Tk() # Creates a first window to choose N the size of the lattice
set_N.title('Select lattice size')
set_N.attributes('-topmost',1) # Brings the window to the front


window = Tk.Tk() # Creates the main window
window.title('Ising 2D')


def quit_all():
   """
   Necessary to quit the GUI properly
   """
   quit_set_N()
   quit_window()
   
def quit_set_N():
    """
    Necessary to quit the GUI properly
    """
    set_N.quit()
    set_N.destroy()

def quit_window():
    """
    Necessary to quit the GUI properly
    """
    window.quit()
    window.destroy()

   
# Creation of a button to quit the whole program.
button_quit = Tk.Button(master=set_N,text='Quit', command=quit_all) 
button_quit.pack(side='bottom',padx=10,pady=10)

# Displays instructions for input
Info = Tk.Label(set_N, text="Enter the size for the lattice \n(integer)")
Info.pack()

# Creation of an input zone
Entry_N = Tk.Entry(set_N, bd=5, width=10,justify='center')
Entry_N.pack()

# Creation of a note
Note = Tk.Text(wrap='word', height=7)
Note.insert('insert', 'Note : To obtain the graphs, it is recommended to use '
                        'a lattice smaller than \nN=64 or computing will' 
                        'become very long and results will be unsatisfying'
                        '\n\nThe recommanded sizes therefore are :\nN=32 for graphs'
                        '\nN=100 for animation')
Note.config(state='disabled')
Note.pack()

# Creation of a button to validate the entry and continue
def _continue():
    """
    Pursues the interaction with user in a second window after N has been selected.
    """
    global N; N=int(Entry_N.get()) # Sets N to the value entered
    quit_set_N() # The input window is quitted
    window.attributes('-topmost', 1) # Brings the main window to the front
       
    # Creation of a 'Quit' button.
    button_quit = Tk.Button(master=window,text='Quit', command=quit_window) 
    button_quit.pack(side='bottom',padx=10,pady=10)
    
    # Creation of a button to obtain the graphs.
    button_graphs = Tk.Button(master=window,text='Obtain Graphs',command=graphs)
    button_graphs.pack(padx=10,pady=10) 
    
    # Creation of a button to launch the animation of the lattice.
    button_animation = Tk.Button(master=window,text='Display animated lattice',command=animated_lattice)
    button_animation.pack(padx=10,pady=10)

button_continue = Tk.Button(master=set_N,text='Continue',command=_continue)
button_continue.pack(side='bottom')

window.mainloop()
set_N.mainloop()
  
