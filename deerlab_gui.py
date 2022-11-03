from tkinter import filedialog
import tkinter
import customtkinter
from inspect import getmembers, signature
import deerlab as dl 
import threading
import matplotlib.pyplot as plt 
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
from itertools import cycle
import time 

running = True 
thread_results = []
customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

from PIL import Image, ImageTk

class App(customtkinter.CTk):

    WIDTH = 1100
    HEIGHT = 700
        
    darker_bckg = '#2a2d2e'
    dark_bckg = '#343638'
    dark_txt = '#afb6be'
    blue = '#1f6aa5'


    dd_modelnames = ['(None) Non-parametric distance distribution']
    dd_models = ['None']
    for model in getmembers(dl):
        model_name, model_obj = model
        if 'dd_' in model_name and model_name!='dd_models': 
            dd_modelnames.append(f'({model_name}) {model_obj.description}') 
            dd_models.append(model_name)
            if 'dd_gauss'==model_name:
                dd_default = model_name
    
    bg_modelnames = []
    bg_models = []
    for model in getmembers(dl):
        model_name, model_obj = model
        if 'bg_' in model_name and model_name!='bg_models': 
            bg_modelnames.append(f'({model_name}) {model_obj.description}') 
            bg_models.append(model_name)
            if 'bg_hom3d'==model_name:
                bg_default = model_name      

    ex_modelnames = []
    ex_models = []
    for model in getmembers(dl):
        model_name, model_obj = model
        if 'ex_' in model_name and model_name!='bg_hom3dex_phase': 
            description = ' '.join(getattr(dl,model_name).__doc__.split('\n')[1].split()[2:])[:-1]
            ex_modelnames.append(f'({model_name}) {description}') 
            ex_models.append(model_name)

    def M_95(self,n=0, top=None, lbl=None):
        # Play GIF (file name = m95.gif) in a 320x320 tkinter window
        # Play GIF concurrently with the loading animation below
        # Close tkinter window after play
        global process_is_alive

        process_is_alive = True
        num_cycles = 8
        count = len(self.frames) * num_cycles
        delay = 8000 // count # make required cycles of animation in around 4 secs
        if n == 0:
            #self.withdraw()
            self.lbl = tkinter.Label(master=self.frame_left, image=self.frames[0], width=160, height=80)
            self.lbl.grid(row=5, column=0)
            process_is_alive = True
            self.lbl.after(delay, self.M_95, n+1, top, self.lbl)
        elif n < count-1:
            self.lbl.config(image=self.frames[n%len(self.frames)])
            self.lbl.after(delay, self.M_95, n+1, top, self.lbl)
        else:
            self.lbl.destroy()
            process_is_alive = False


    def plot_distribution(self):
            
        if hasattr(self,'frame_distribution'):
            self.frame_distribution.destroy()
        
        self.frame_distribution = customtkinter.CTkFrame(master=self.frame_right,width=400,height=250,fg_color='#ffffff')
        self.frame_distribution.grid(row=0, column=1)
        
        # the figure that will contain the plot
        px2in = 0.0104
        figwidth = self.frame_dataplot.cget('width')
        figheight = self.frame_dataplot.cget('height')
        fig = Figure(figsize = (px2in*figwidth, px2in*figheight),
                    facecolor=App.darker_bckg)

        # adding the subplot
        ax = fig.add_subplot(111)
        fig.set_tight_layout(True)
        # Change the background color
        ax.set_facecolor(App.darker_bckg)

        if hasattr(self,'results'): 

            # Plot the data
            ax.plot(self.results['r'],self.results['P'],'-',color=App.blue, markersize=5)
            ax.fill_between(self.results['r'],*self.results['PUncert'].ci(95).T,linewidth=0,color=App.blue, alpha=0.5)
            ax.set_ylabel('P(r) [1/nm]')
            ax.set_xlabel('r [nm]')

            # Make the top and right spines invisible
            ax.spines[['top', 'right']].set_visible(False)

            # Set color of axes
            for axis in ['bottom','top','right','left']:
                ax.spines[axis].set_color(App.dark_txt)
                ax.spines[axis].set_linewidth(1.5)
            ax.xaxis.label.set_color(App.dark_txt)
            ax.yaxis.label.set_color(App.dark_txt)
            ax.tick_params(colors=App.dark_txt, which='both') 
        else: 
            ax.text(0.4,0.5,'No results',color=App.dark_txt)
            ax.axis('off')
        # Creating the Tkinter canvas containing the Matplotlib figure
        self.data_plot = FigureCanvasTkAgg(fig, master=self.frame_distribution)  
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack()


    def run_analysis(self):

        global process_is_alive

        ex_model = getattr(dl,[ex_model for ex_model in App.ex_models if ex_model in self.Exmodel_menu.get()][0])     
        bg_model = getattr(dl,[bg_model for bg_model in App.bg_models if bg_model in self.Bmodel_menu.get()][0])        
        dd_model = [dd_model for dd_model in App.dd_models if dd_model in self.Pmodel_menu.get()][0]   
        if dd_model=='None':
            dd_model = None
        else: 
            dd_model = getattr(dl,dd_model)

        # Construct distance vector
        rmin = float(self.rmin_entry.get())            
        rmax = float(self.rmax_entry.get())            
        dr = float(self.dr_entry.get())            
        r = np.arange(rmin,rmax,dr)

        t = self.data['t']
        t = t - t[0] + float(self.deadtime_entry.get())
        self.data['t'] = t

        # Construct experiment model 
        taus = [float(getattr(self,f"delay{n+1}_entry").get()) for n in range(self.Ndelays)]
        pathways = [n+1 for n in range(self.Npathways) if bool(getattr(self,f'pathway{n+1}_switch').get())]
        experiment = ex_model(*taus,pathways=pathways)

        Vmodel = dl.dipolarmodel(self.data['t'],r,Pmodel=dd_model, Bmodel=bg_model, experiment=experiment)

        cancel_id = None
        loading_label = customtkinter.CTkLabel(master=self.frame_left, width=16, height=8, bg_color=App.darker_bckg)
        loading_label.grid(row=5, column=0)
        def start_loading():
            global running
            global thread_results
            if running:  # Animation not started?
                loading_label.configure(image=next(self.frames))
                app.after(10, start_loading) # call this function every 100ms
            else:  
                loading_label.destroy()       
                results = thread_results[0]
                if not hasattr(results,'P'):
                    P = results.evaluate(dd_model,r)
                    PUncert = results.propagate(dd_model,r,lb=np.zeros_like(r)) 
                else: 
                    P = results.P 
                    PUncert = results.PUncert

                self.results = {'r':r,'P':P,'PUncert':PUncert,'t':self.data['t'], 'model':results.model}
                self.plot_distribution()
                self.plot_data()

        def loadingAnimation():
            global running
            global thread_results
            running = True
            print('Running')
            results = dl.fit(Vmodel,self.data['data'])
            print('Finished')
            running = False
            thread_results.append(results) 
 
        threading.Thread(target=loadingAnimation).start()
        start_loading()

        return

    def plot_data(self):

            
        if hasattr(self,'frame_dataplot'):
            self.frame_dataplot.destroy()
        
        self.frame_dataplot = customtkinter.CTkFrame(master=self.frame_right,width=400,height=250,fg_color='#ffffff')
        self.frame_dataplot.grid(row=0, column=0, sticky="")
        
        # the figure that will contain the plot
        px2in = 0.0104
        figwidth = self.frame_dataplot.cget('width')
        figheight = self.frame_dataplot.cget('height')
        fig = Figure(figsize = (px2in*figwidth, px2in*figheight),
                    facecolor=App.darker_bckg)

        # adding the subplot
        ax = fig.add_subplot(111)
        fig.set_tight_layout(True)
        # Change the background color
        ax.set_facecolor(App.darker_bckg)

        if hasattr(self,'data'): 

            # Plot the data
            ax.plot(self.data['t'],self.data['data'],'.',color=App.dark_txt, markersize=5)

            if hasattr(self,'results'): 
                ax.plot(self.results['t'],self.results['model'],'-',color=App.blue, linewidth=3)

            ax.set_ylabel('V(t) [arb.u.]')
            ax.set_xlabel('t [us]')

            # Make the top and right spines invisible
            ax.spines[['top', 'right']].set_visible(False)

            # Set color of axes
            for axis in ['bottom','top','right','left']:
                ax.spines[axis].set_color(App.dark_txt)
                ax.spines[axis].set_linewidth(1.5)
            ax.xaxis.label.set_color(App.dark_txt)
            ax.yaxis.label.set_color(App.dark_txt)
            ax.tick_params(colors=App.dark_txt, which='both') 
        else: 
            ax.text(0.4,0.5,'No data',color=App.dark_txt)
            ax.axis('off')
        
        # Creating the Tkinter canvas containing the Matplotlib figure
        self.data_plot = FigureCanvasTkAgg(fig, master=self.frame_dataplot)  
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack()


    def load_file(self):
        file = filedialog.askopenfilename()
        t,Vexp = dl.deerload(file)
        Vexp = dl.correctphase(Vexp,offset=True) 
        Vexp = Vexp/max(Vexp)
        t = t - t[0]
        self.data = {'t':t, 'data':Vexp}
        self.plot_data()
        self.run_button.configure(state='normal')


    def setup_pulsedelays(self,ex_modelname):
        for n in range(5):
            try:
                getattr(self,f"delay{n+1}_label").destroy()
                getattr(self,f"delay{n+1}_entry").destroy()
            except: pass
        ex_model = [ex_model for ex_model in App.ex_models if ex_model in ex_modelname ][0]        
        sig = signature(getattr(dl,ex_model))
        Ndelays = len([param for param in sig.parameters.values() if 'tau' in param.name])
        self.Ndelays = Ndelays
        self.deadtime_label = customtkinter.CTkLabel(master=self.frame_pulsedelays, text=f"Dead time",width=50)
        self.deadtime_label.grid(row=0, column=0,padx=5,sticky='we')
        self.deadtime_entry = customtkinter.CTkEntry(master=self.frame_pulsedelays, placeholder_text="us",width=50)
        self.deadtime_entry.grid(row=1, column=0,padx=5,sticky='we')
        ncol = 1
        for n in range(Ndelays):
            setattr(self,f"delay{n+1}_label", customtkinter.CTkLabel(master=self.frame_pulsedelays, text=f"tau{n+1}",width=50) )
            getattr(self,f"delay{n+1}_label").grid(row=0, column=ncol,padx=5,sticky='we')
            setattr(self,f"delay{n+1}_entry", customtkinter.CTkEntry(master=self.frame_pulsedelays, placeholder_text="us",width=50) )
            getattr(self,f"delay{n+1}_entry").grid(row=1, column=ncol,padx=5,sticky='we')
            ncol += 1


    def setup_pathways(self,ex_modelname):
        ex_model = [ex_model for ex_model in App.ex_models if ex_model in ex_modelname ][0]        
        sig = signature(getattr(dl,ex_model))
        Ndelays = len([param for param in sig.parameters.values() if 'tau' in param.name])
        Npathways = getattr(dl,ex_model)(*[0]*Ndelays).npathways
        self.Npathways = Npathways
        self.Ndelays = Ndelays

        for n in range(10):
            try:
                getattr(self,f"pathway{n+1}_switch").destroy()
            except: pass

        nrow,ncol = 0,0
        for n in range(Npathways):
            setattr(self,f"pathway{n+1}_switch", customtkinter.CTkSwitch(master=self.frame_pathways, text=f"#{n+1}", command=None, onvalue=True, offvalue=False))
            getattr(self,f"pathway{n+1}_switch").grid(row=nrow, column=ncol, padx=10)
            ncol+=1
            if ncol>3:
                nrow = 1
                ncol = 0
        self.pathway1_switch.select()

    def change_experiment(self,ex_model):
        self.setup_pulsedelays(ex_model)
        self.setup_pathways(ex_model)

    def __init__(self):
        super().__init__()

        self.title("CustomTkinter complex_example.py")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # call .on_closing() when app gets closed


        frames = []
        im = Image.open(r"D:\lufa\projects\DeerLab\testing\GIF4.gif")
        for frame in range(im.n_frames):
            image = im.seek(frame)
            # For each pixel in the image
            for i in range(image.size[0]):
                for j in range(image.size[1]):
                    # If the pixel is white
                    if pixels[i, j] == (255, 255, 255, 255):
                        # Make it transparent
                        pixels[i, j] = (255, 255, 255, 0)
            frames.append(ImageTk.PhotoImage(im.convert('RGBA').resize((160, 80))))
        frames_cycle=cycle(frames)
        self.frames = frames_cycle

        # ============ create two frames ============

        # configure grid layout (2x1)
        self.grid_columnconfigure((0,1), weight=1)
        self.grid_rowconfigure((0,1), weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self,
                                                 width=180,
                                                 corner_radius=0)
        self.frame_left.grid(row=0, column=0)

        self.frame_right = customtkinter.CTkFrame(master=self,width=800)
        self.frame_right.grid(row=0, column=1, sticky="nswe", padx=20, pady=20)

        # ============ frame_left ============

        # configure grid layout (1x11)
        self.frame_left.grid_rowconfigure(0, minsize=10)   # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(5, weight=1)     # empty row as spacing
        self.frame_left.grid_rowconfigure(8, minsize=20)    # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(11, minsize=10)  # empty row with minsize as spacing

        self.label_1 = customtkinter.CTkLabel(master=self.frame_left,
                                              text="CustomTkinter",
                                              text_font=("Roboto Medium", -16))  # font name and size in px
        self.label_1.grid(row=1, column=0, pady=10, padx=10)

        self.button_1 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Load file",
                                                command=self.load_file)
        self.button_1.grid(row=2, column=0, pady=10, padx=20)

        self.button_2 = customtkinter.CTkButton(master=self.frame_left,
                                                state = 'disabled',
                                                text="Generate script",
                                                command=self.button_event)
        self.button_2.grid(row=3, column=0, pady=10, padx=20)

        self.run_button = customtkinter.CTkButton(master=self.frame_left,
                                                state = 'disabled',
                                                text="Run analysis",
                                                command=self.run_analysis)
        self.run_button.grid(row=4, column=0, pady=10, padx=20)

        self.label_mode = customtkinter.CTkLabel(master=self.frame_left, text="Appearance Mode:")
        self.label_mode.grid(row=9, column=0, pady=0, padx=20, sticky="w")

        self.optionmenu_1 = customtkinter.CTkOptionMenu(master=self.frame_left,
                                                        values=["Light", "Dark", "System"],
                                                        command=self.change_appearance_mode)
        self.optionmenu_1.grid(row=10, column=0, pady=10, padx=20, sticky="w")

        # ============ frame_right ============

        # configure grid layout (2x2)
        self.frame_right.grid_rowconfigure((0, 1), weight=1)
        self.frame_right.grid_columnconfigure((0, 1), weight=1)


        self.frame_modelling = customtkinter.CTkFrame(master=self.frame_right,width=400,fg_color=App.dark_bckg)
        self.frame_modelling.grid(row=1, column=0, sticky="", pady=15)

        self.frame_analysis = customtkinter.CTkFrame(master=self.frame_right,width=400,fg_color=App.darker_bckg)
        self.frame_analysis.grid(row=1, column=1, sticky="")

        self.plot_data()
        self.plot_distribution()
        # ============ frame_modelling ============
        self.frame_modelling.grid_rowconfigure(9, weight=1)
        self.frame_modelling.grid_columnconfigure(0, weight=1)
        self.frame_modelling.grid_columnconfigure(1, weight=1000)
        
        self.Modelling_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Modelling", text_font='Helvetica 14 bold')
        self.Modelling_label.grid(row=0, column=0, columnspan=2, pady=10, sticky="we")

        self.Pmodel_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Distribution:",width=50)
        self.Pmodel_label.grid(row=1, column=0, pady=5, padx=15, sticky="we")
        self.Pmodel_menu = customtkinter.CTkOptionMenu(self.frame_modelling,width=250, values=App.dd_modelnames,dynamic_resizing=False)
        self.Pmodel_menu.set(App.dd_modelnames[0])
        self.Pmodel_menu.grid(row=1, column=1, pady=5, padx=10, sticky="we")

        self.Bmodel_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Background:",width=50)
        self.Bmodel_label.grid(row=2, column=0, pady=5, padx=15, sticky="we")
        self.Bmodel_menu = customtkinter.CTkOptionMenu(self.frame_modelling,width=250, values=App.bg_modelnames,dynamic_resizing=False)
        self.Bmodel_menu.set(App.bg_modelnames[1])
        self.Bmodel_menu.grid(row=2, column=1, pady=5, padx=10, sticky="we")

        self.Exmodel_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Experiment:",width=50)
        self.Exmodel_label.grid(row=3, column=0, pady=5, padx=15, sticky="we")
        self.Exmodel_menu = customtkinter.CTkOptionMenu(self.frame_modelling,width=250, values=App.ex_modelnames,dynamic_resizing=False,command=self.change_experiment)
        self.Exmodel_menu.set(App.ex_modelnames[1])
        self.Exmodel_menu.grid(row=3, column=1, pady=5, padx=10, sticky="we")


        self.frame_pulsedelays = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_pulsedelays.grid(row=5, column=0, columnspan=2, pady=4, sticky="ns")


        self.setup_pulsedelays('ex_4pdeer') 

        self.pathways_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Dipolar pathways")
        self.pathways_label.grid(row=6, column=0, columnspan=2, pady=2, sticky="we")

        self.frame_pathways = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_pathways.grid(row=7, column=0, columnspan=2)

        self.setup_pathways('ex_4pdeer') 


        self.frame_distances = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_distances.grid(row=9, column=0, columnspan=2)


        self.distances_label = customtkinter.CTkLabel(master=self.frame_distances, text="Distance range",width=70)
        self.distances_label.grid(row=1, column=0, rowspan=2,  pady=2, sticky="we")

        self.rmin_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"min.",width=70) 
        self.rmin_label.grid(row=0, column=1,padx=5,sticky='we')
        self.rmin_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=70)
        self.rmin_entry.grid(row=1, column=1,padx=5,sticky='we')
        self.rmax_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"max.",width=70) 
        self.rmax_label.grid(row=0, column=2,padx=5,sticky='we')
        self.rmax_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=70)
        self.rmax_entry.grid(row=1, column=2,padx=5,sticky='we')
        self.dr_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"resolution",width=70) 
        self.dr_label.grid(row=0, column=3,padx=5,sticky='wes')
        self.dr_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=70)
        self.dr_entry.grid(row=1, column=3,padx=5,sticky='we')

        self.rmin_entry.insert(0,'1.5')
        self.rmax_entry.insert(0,'8')
        self.dr_entry.insert(0,'0.05')
    def button_event(self):
        print("Button pressed")

    def change_appearance_mode(self, new_appearance_mode):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def on_closing(self, event=0):
        self.destroy()


if __name__ == "__main__":
    app = App()
    app.mainloop()