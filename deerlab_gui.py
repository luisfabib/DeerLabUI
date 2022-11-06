from tkinter import filedialog
import tkinter
import customtkinter
from inspect import getmembers, signature
import deerlab as dl 
import threading
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
from itertools import cycle
from PIL import Image, ImageTk
from webbrowser import open_new
import os 

# Path to current file
PATH = os.path.dirname(os.path.realpath(__file__))

# Global variables for threading
THREAD_ACTIVE = True 
THREAD_RETURNS = []

# App theme
customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("green")  # Themes: "blue" (standard), "green", "dark-blue"


class App(customtkinter.CTk):

    WIDTH = 1350
    HEIGHT = 800
        
    darker_bckg = '#232937'
    dark_bckg = '#393e4b'
    dark_txt = '#afb6be'

    darker_bckg = '#2a2d2e'
    dark_bckg = '#343638'
    dark_txt = '#afb6be'
    blue = '#1f6aa5'
    magenta = '#df4577'
    orange = '#d7a064'
    light_green = '#27e8a7'
    green = '#11b384'
    dark_green = '#43675b'

    switchcolor = green

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


    #==============================================================================================================
    def plot_distribution(self):
            
        if hasattr(self,'frame_distribution'):
            self.frame_distribution.destroy()
        
        self.frame_distribution = customtkinter.CTkFrame(master=self.frame_right,width=650,height=250,fg_color=App.darker_bckg)
        self.frame_distribution.grid(row=0, column=1, padx=5)
        
        # the figure that will contain the plot
        px2in = 0.0104
        figwidth = self.frame_distribution.cget('width')
        figheight = self.frame_distribution.cget('height')
        fig = Figure(figsize = (px2in*figwidth, px2in*figheight),
                    facecolor=App.darker_bckg)

        # adding the subplot
        ax = fig.add_subplot(111)
        fig.set_tight_layout(True)
        # Change the background color
        ax.set_facecolor(App.darker_bckg)

        if hasattr(self,'results'): 

            # Plot the data
            ax.fill_between(self.results['r'],*self.results['PUncert'].ci(95).T,linewidth=0,color=App.light_green, alpha=0.4)
            ax.fill_between(self.results['r'],*self.results['PUncert'].ci(50).T,linewidth=0,color=App.light_green, alpha=0.4)
            ax.plot(self.results['r'],self.results['P'],'-', linewidth=2,color=App.dark_green)
            ax.set_xlabel('r [nm]')

            # Make the top and right spines invisible
            ax.spines[['top', 'right','left']].set_visible(False)

            # Set color of axes
            for axis in ['bottom']:
                ax.spines[axis].set_color(App.dark_txt)
                ax.spines[axis].set_linewidth(1.5)
            ax.xaxis.label.set_color(App.dark_txt)
            ax.tick_params(left = False, right = False , labelleft = False)
            ax.tick_params(colors=App.dark_txt, which='both') 
            ax.tick_params(colors=App.dark_txt, which='both') 
            ax.autoscale(enable=True, axis='both', tight=True)
        else: 
            ax.text(0.4,0.5,'No results',color=App.dark_txt)
            ax.axis('off')
        # Creating the Tkinter canvas containing the Matplotlib figure
        self.data_plot = FigureCanvasTkAgg(fig, master=self.frame_distribution)  
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack()
    #==============================================================================================================



    #==============================================================================================================
    def run_analysis(self):

        global THREAD_ACTIVE, THREAD_RETURNS

        # Disable main menu buttons during execution
        self.run_button.configure(text='Running...', state='disabled')
        self.load_button.configure(state='disabled')

        # Get the selection of models 
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

        # Adjust the time vector by the specified deadtime
        t = self.data['t']
        t = t - t[0] + float(self.deadtime_entry.get())
        self.data['t'] = t

        # Construct the experiment model 
        taus = [float(getattr(self,f"delay{n+1}_entry").get()) for n in range(self.Ndelays)]
        pathways = [n+1 for n in range(self.Npathways) if bool(getattr(self,f'pathway{n+1}_switch').get())]
        experiment = ex_model(*taus,pathways=pathways)

        # Construct the dipolar model
        Vmodel = dl.dipolarmodel(self.data['t'],r,Pmodel=dd_model, Bmodel=bg_model, experiment=experiment)

        # Get options for the analysis 
        regparam = self.regparam_menu.get().lower()
        if self.compactness_switch.get():
            compactness = dl.dipolarpenalty(dd_model,r,'compactness')
        else:
            compactness = None
        if self.bootstrap_switch.get():
            bootstrap = int(self.bootstrap_entry.get())
        else:
            bootstrap = 0
        # Construct a canvas to put the runtime animation
        animation_canvas = tkinter.Canvas(master=self.frame_left, width=180, height=60, bg=App.darker_bckg, highlightthickness=0)
        animation_canvas.grid(row=5, column=0)

        #--------------------------------------------------------------------------------------------------
        def animated_wait_for_results():

            global THREAD_ACTIVE, THREAD_RETURNS

            # Check if the analysis thread is still active 
            if THREAD_ACTIVE:
                # If still running, keep on updating the animation
                animation_canvas.create_image(0,0,image=next(self.frames),anchor='nw')
                app.after(50, animated_wait_for_results) # Animation update every 50ms
            else:  
                # If finished, stop the animation and retrieve the output of the analysis 
                animation_canvas.destroy()       
                results = THREAD_RETURNS[0]
                
                # Evaluate the distance distribution estimate and uncertainty
                if not hasattr(results,'P'):
                    P = results.evaluate(dd_model,r)
                    PUncert = results.propagate(dd_model,r,lb=np.zeros_like(r)) 
                else: 
                    P = results.P 
                    PUncert = results.PUncert

                # Clear output console
                self.textbox.delete('1.0', tkinter.END)
                # Print results summary table
                self.textbox.insert(tkinter.INSERT,results._summary)

                # Update the plots with the analysis results
                self.results = {'r':r,'P':P,'PUncert':PUncert,'t':self.data['t'], 'model':results.model}
                self.plot_distribution()
                self.plot_data()

                # Reactivate the main menu buttons
                self.run_button.configure(text='Run analysis', state='normal')
                self.load_button.configure(state='normal')
                self.script_button.configure(state='normal')
                return # Finish the analysis and return to mainloop()
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        def threaded_analysis():
            
            global THREAD_ACTIVE, THREAD_RETURNS
            
            # Set thread status as active
            THREAD_ACTIVE = True
            
            # Fit the dipolar model to the data
            results = dl.fit(Vmodel,self.data['data'], 
                        regparam=regparam, penalties=compactness, bootstrap=bootstrap)
            
            # Pack results to be extracted outside the thread
            THREAD_RETURNS.append(results) 

            # Set status of the thread to inactive before killing it
            THREAD_ACTIVE = False
        #--------------------------------------------------------------------------------------------------
    
        # Prepare container for thread outputs
        THREAD_RETURNS = []
        # Start analysis on a separate thread
        threading.Thread(target=threaded_analysis).start()
        # Run an animation while waiting for the results to finish 
        animated_wait_for_results()
    #==============================================================================================================


    #==============================================================================================================
    def plot_data(self):

            
        if hasattr(self,'frame_dataplot'):
            self.frame_dataplot.destroy()
        
        self.frame_dataplot = customtkinter.CTkFrame(master=self.frame_right,width=400,height=250,fg_color=App.darker_bckg)
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
            ax.plot(self.data['t'],self.data['data'],'.',color=App.dark_green, markersize=5)

            if hasattr(self,'results'): 
                ax.plot(self.results['t'],self.results['model'],'-',color=App.green, linewidth=2.5)

            ax.set_ylabel('V(t) [arb.u.]')
            ax.set_xlabel('t [μs]')

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
    #==============================================================================================================


    #==============================================================================================================
    def load_file(self):

        # Us the OS dialog window to select a file
        file = filedialog.askopenfilename()
        if file=='': return
        self.filepath = file 

        # Load the file with DeerLab
        t,Vexp = dl.deerload(file)

        # Phase correction
        Vexp = dl.correctphase(Vexp,offset=True) 
        Vexp = Vexp/max(Vexp)

        # Adjust time axis
        t = t - t[0]

        # If there are already results from a previous dataset, delete them
        if hasattr(self,'results'):
            delattr(self,'results')
 
        # Store experimental data into the app
        self.data = {'t':t, 'data':Vexp}

        # Update the data and results displays
        self.plot_data()
        self.plot_distribution()

        # Enable the dependent buttons 
        self.run_button.configure(state='normal')
        self.autodistances_button.configure(state='normal')

    #==============================================================================================================


    #==============================================================================================================
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
    #==============================================================================================================


    #==============================================================================================================
    def setup_pathways(self,ex_modelname):
        
        Npathways_dict = {
            'ex_3pdeer': 2,
            'ex_4pdeer': 4,
            'ex_fwd5pdeer': 8,
            'ex_rev5pdeer': 8,
            'ex_dqc': 8,
            'ex_sifter': 3, 
            'ex_ridme': 4
        }

        ex_model = [ex_model for ex_model in App.ex_models if ex_model in ex_modelname ][0]        
        sig = signature(getattr(dl,ex_model))
        Ndelays = len([param for param in sig.parameters.values() if 'tau' in param.name])
        Npathways = Npathways_dict[ex_model]
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
    #==============================================================================================================


    #==============================================================================================================
    def change_experiment(self,ex_model):
        self.setup_pulsedelays(ex_model)
        self.setup_pathways(ex_model)
    #==============================================================================================================



    #==============================================================================================================
    def __init__(self):
        super().__init__()

        self.title("DeerLab UI - Dipolar EPR spectroscopy")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")
        self.configure(bg='#151921')
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # call .on_closing() when app gets closed
                
        # Load the animation GIF into the UI's memory
        frames = []
        im = Image.open(PATH + "\graphics\loading.gif")
        for frame in range(im.n_frames):
            # Get current frame
            im.seek(frame)
            # Process frame
            frames.append(ImageTk.PhotoImage(im.convert('RGBA').resize((180, 60))))
        # Construct looped generator for the GIF frames
        frames_cycle=cycle(frames)
        self.frames = frames_cycle

        # ============ create two frames ============

        # configure grid layout (2x1)
        self.grid_columnconfigure((0,1), weight=1)
        self.grid_rowconfigure((0), weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self,
                                                 width=230,
                                                 fg_color=App.darker_bckg,
                                                 corner_radius=5)
        self.frame_left.grid(row=0, column=0, padx=(20,10),sticky="",)

        self.frame_right = customtkinter.CTkFrame(master=self,width=750,
                                                 fg_color=App.darker_bckg,)
        self.frame_right.grid(row=0, column=1, sticky="ns", padx=(0,0), pady=(10,10))

        # ============ frame_left ============

        # configure grid layout (1x11)
        self.frame_left.grid_rowconfigure(0, minsize=10)   # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(5, weight=1)     # empty row as spacing
        self.frame_left.grid_rowconfigure(8, minsize=20)    # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(11, minsize=10)  # empty row with minsize as spacing

        self.label_1 = customtkinter.CTkLabel(master=self.frame_left,
                                              text="Main Menu",
                                              text_font=("Roboto Medium", -16))  # font name and size in px
        self.label_1.grid(row=1, column=0, pady=10, padx=10)

        self.load_button_image = self.load_image("/graphics/folder.png", 30)
        self.load_button = customtkinter.CTkButton(master=self.frame_left,
                                                text="Load dataset",
                                                height=40,
                                                image=self.load_button_image,
                                                compound="right",
                                                fg_color=App.magenta,
                                                hover_color='#7d3634',
                                                command=self.load_file)
        self.load_button.grid(row=2, column=0, pady=10, padx=20)

        self.script_button_image = self.load_image("/graphics/script.png", 40)
        self.script_button = customtkinter.CTkButton(master=self.frame_left,
                                                state = 'disabled',
                                                text="Script",
                                                height=40,
                                                fg_color='#656565',
                                                image=self.script_button_image,
                                                compound="right",
                                                command=self.generate_script)
        self.script_button.grid(row=11, column=0, pady=10, padx=20)


        self.report_button_image = self.load_image("/graphics/report.png", 30)
        self.report_button = customtkinter.CTkButton(master=self.frame_left,
                                                state = 'disabled',
                                                text="Report",
                                                height=40,
                                                fg_color='#656565',
                                                image=self.report_button_image,
                                                compound="right",
                                                command=None)
        self.report_button.grid(row=12, column=0, pady=10, padx=20)

        self.run_button_image = self.load_image("/graphics/run.png", 30)
        self.run_button = customtkinter.CTkButton(master=self.frame_left,
                                                state = 'disabled',
                                                text="Run analysis",
                                                height=40,
                                                fg_color=App.green,
                                                image=self.run_button_image,
                                                compound="right",
                                                command=self.run_analysis)
        self.run_button.grid(row=3, column=0, pady=10, padx=20)

        #self.label_mode = customtkinter.CTkLabel(master=self.frame_left, text="Appearance Mode:")
        #self.label_mode.grid(row=9, column=0, pady=0, padx=20, sticky="w")
        #self.optionmenu_1 = customtkinter.CTkOptionMenu(master=self.frame_left,values=["Light", "Dark", "System"],command=self.change_appearance_mode)
        #self.optionmenu_1.grid(row=10, column=0, pady=10, padx=20, sticky="w")

        # ============ frame_right ============

        # configure grid layout (2x2)
        self.frame_right.grid_rowconfigure((0, 1), weight=1)
        self.frame_right.grid_columnconfigure((0, 1), weight=1)

        self.frame_modelling = customtkinter.CTkFrame(master=self.frame_right,width=400,fg_color=App.dark_bckg)
        self.frame_modelling.grid(row=1, column=0, sticky="nswe", pady=(0,10),  padx=10)

        self.plot_data()
        self.plot_distribution()

        # ============ frame_modelling ============
        self.frame_modelling.grid_rowconfigure(20, weight=1)
        self.frame_modelling.grid_columnconfigure(0, weight=1)
        self.frame_modelling.grid_columnconfigure(1, weight=1000)
        
        self.Modelling_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Modelling", text_font='Helvetica 13 bold')
        self.Modelling_label.grid(row=0, column=0, pady=(10,0), sticky="we")

        self.frame_models = customtkinter.CTkFrame(master=self.frame_modelling,width=400,fg_color=App.dark_bckg)
        self.frame_models.grid(row=1, column=0, sticky="we", pady=0,  padx=(15,0))

        self.Exmodel_label = customtkinter.CTkLabel(master=self.frame_models, text="Experiment",width=50,text_font='Helvetica 10 bold')
        self.Exmodel_label.grid(row=2, column=0, pady=5, sticky="we")
        self.Exmodel_menu = customtkinter.CTkOptionMenu(self.frame_models,width=280, values=App.ex_modelnames,dynamic_resizing=False, command=self.change_experiment)
        self.Exmodel_menu.set(App.ex_modelnames[1])
        self.Exmodel_menu.grid(row=2, column=1, pady=(5,0), padx=5, sticky="we")

        self.frame_pulsedelays = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_pulsedelays.grid(row=3, column=0, columnspan=2, pady=(0,5), sticky="ns")
        self.setup_pulsedelays('ex_4pdeer') 

        self.frame_pathways_label = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_pathways_label.grid(row=4, column=0, columnspan=2, pady=(0,5), sticky="ns")
        self.frame_modelling.grid_columnconfigure(0, weight=1000)
        self.pathways_label = customtkinter.CTkLabel(master=self.frame_pathways_label, text="Dipolar pathways")
        self.pathways_label.grid(row=0, column=0, pady=2, sticky="we")
        self.pathways_label_help = customtkinter.CTkButton(master=self.frame_pathways_label, width=7, height=7, text="?", command=self.openbrowser, fg_color='#656565')
        self.pathways_label_help.grid(row=0, column=2, pady=2, sticky="w")
        self.frame_pathways = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_pathways.grid(row=5, column=0, columnspan=2)
        self.setup_pathways('ex_4pdeer') 

        self.frame_Pmodels = customtkinter.CTkFrame(master=self.frame_modelling,width=400,fg_color=App.dark_bckg)
        self.frame_Pmodels.grid(row=6, column=0, sticky="we", pady=(10,0),  padx=15)

        self.Pmodel_label = customtkinter.CTkLabel(master=self.frame_Pmodels, text="Distribution",width=50,text_font='Helvetica 10 bold')
        self.Pmodel_label.grid(row=0, column=0, pady=5, sticky="we")
        self.Pmodel_menu = customtkinter.CTkOptionMenu(self.frame_Pmodels,width=280, values=App.dd_modelnames,dynamic_resizing=False, )
        self.Pmodel_menu.set(App.dd_modelnames[0])
        self.Pmodel_menu.grid(row=0, column=1, pady=5, padx=5, sticky="we")

        self.frame_distances = customtkinter.CTkFrame(master=self.frame_modelling,fg_color=App.dark_bckg)
        self.frame_distances.grid(row=7, column=0, columnspan=2)
        self.distances_label = customtkinter.CTkLabel(master=self.frame_distances, text="Distance range",width=70)
        self.distances_label.grid(row=1, column=0, rowspan=2,  pady=2, sticky="we")
        self.rmin_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"min.",width=50) 
        self.rmin_label.grid(row=0, column=1,padx=5,sticky='we')
        self.rmin_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=50)
        self.rmin_entry.grid(row=1, column=1,padx=5,sticky='we')
        self.rmax_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"max.",width=50) 
        self.rmax_label.grid(row=0, column=2,padx=5,sticky='we')
        self.rmax_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=50)
        self.rmax_entry.grid(row=1, column=2,padx=5,sticky='we')
        self.dr_label = customtkinter.CTkLabel(master=self.frame_distances, text=f"resolution",width=50) 
        self.dr_label.grid(row=0, column=3,padx=5,sticky='wes')
        self.dr_entry = customtkinter.CTkEntry(master=self.frame_distances, placeholder_text="nm",width=50)
        self.dr_entry.grid(row=1, column=3,padx=5,sticky='we')
        self.autodistances_button = customtkinter.CTkButton(master=self.frame_distances,state = 'disableds', width=50, text="auto", command=self.automatic_distances)
        self.autodistances_button.grid(row=1, column=4, rowspan=2, padx=5, sticky="we")
        self.rmin_entry.insert(0,'1.5')
        self.rmax_entry.insert(0,'8')
        self.dr_entry.insert(0,'0.05')

        self.frame_Bmodels = customtkinter.CTkFrame(master=self.frame_modelling,width=400,fg_color=App.dark_bckg)
        self.frame_Bmodels.grid(row=8, column=0, sticky="we", pady=(10,0),  padx=15)

        self.Bmodel_label = customtkinter.CTkLabel(master=self.frame_Bmodels, text="Background",width=50,text_font='Helvetica 10 bold')
        self.Bmodel_label.grid(row=1, column=0, pady=5, sticky="we")
        self.Bmodel_menu = customtkinter.CTkOptionMenu(self.frame_Bmodels,width=280, values=App.bg_modelnames, dynamic_resizing=False)
        self.Bmodel_menu.set(App.bg_modelnames[1])
        self.Bmodel_menu.grid(row=1, column=1, pady=5, padx=5, sticky="we")

        self.Analysis_label = customtkinter.CTkLabel(master=self.frame_modelling, text="Analysis", text_font='Helvetica 13 bold')
        self.Analysis_label.grid(row=10, column=0, columnspan=2, pady=5, sticky="")

        self.regparam_frame = customtkinter.CTkFrame(master=self.frame_modelling,width=300,fg_color=App.dark_bckg)
        self.regparam_frame.grid(row=11, column=0, pady=0, sticky="")

        self.regparam_label = customtkinter.CTkLabel(master=self.regparam_frame, text="Smoothness regularization:",width=80)
        self.regparam_label.grid(row=0, column=0, pady=0, sticky="we")
        self.regparam_menu = customtkinter.CTkOptionMenu(self.regparam_frame,width=80, values=['AIC','BIC','cAIC','GCV','srGCV','LR','LC'],
                        dynamic_resizing=False, button_color=App.green,fg_color=App.dark_bckg, dropdown_hover_color=App.green)
        self.regparam_menu.set('AIC')
        self.regparam_menu.grid(row=0, column=1, pady=0, padx=5, sticky="we")

        self.compactness_switch = customtkinter.CTkSwitch(master=self.frame_modelling, text="Compactness regularization", onvalue=True, offvalue=False)
        self.compactness_switch.grid(row=12, column=0, columnspan=2, padx=(0,35), pady=5, sticky="")

        self.bootstrap_frame = customtkinter.CTkFrame(master=self.frame_modelling,width=300,fg_color=App.dark_bckg)
        self.bootstrap_frame.grid(row=13, column=0, pady=0, sticky="")
        self.bootstrap_switch = customtkinter.CTkSwitch(master=self.bootstrap_frame, command=self.bootstrap_switch_samples, text="Bootstrapping", onvalue=True, offvalue=False)
        self.bootstrap_switch.grid(row=0, column=0, pady=0, sticky="")
        self.bootstrap_entry = customtkinter.CTkEntry(master=self.bootstrap_frame, width=50, state='disabled')
        self.bootstrap_entry.grid(row=0, column=1,  padx=5, pady=0, sticky="")
        self.bootstrap_label = customtkinter.CTkLabel(master=self.bootstrap_frame, text=f"Samples",width=50) 
        self.bootstrap_label.grid(row=0, column=2,pady=0,sticky='')
        


        self.frame_results = customtkinter.CTkFrame(master=self.frame_right,width=650,fg_color=App.dark_bckg)
        self.frame_results.grid(row=1, column=1, sticky="news", pady=(0,10),  padx=10)
        self.frame_results.grid_columnconfigure(0, weight=1)
        self.frame_results.grid_rowconfigure((1), weight=10)


        self.Results_label = customtkinter.CTkLabel(master=self.frame_results, text="Results", text_font='Helvetica 13 bold')
        self.Results_label.grid(row=0, column=0, pady=(10,0), sticky="n")

        self.frame_textbox = customtkinter.CTkFrame(master=self.frame_results,width=650,fg_color=App.dark_bckg)
        self.frame_textbox.grid(row=1, column=0, sticky="news")
        self.frame_textbox.grid_rowconfigure(0, weight=10)
        self.frame_textbox.grid_columnconfigure(0, weight=1)

        # create scrollable textbox
        self.textbox = tkinter.Text(self.frame_textbox, highlightthickness=0,bg=App.dark_bckg,font='Consolas 10', relief='flat', fg=App.dark_txt, wrap='none')
        self.textbox.grid(row=0, column=0, padx=10, sticky="news")
        self.textbox_scrollbar = customtkinter.CTkScrollbar(self.frame_textbox, command=self.textbox.yview)
        self.textbox_scrollbar.grid(row=0, column=1, sticky="ns")
        self.textbox_scrollbarx = customtkinter.CTkScrollbar(self.frame_textbox, command=self.textbox.xview, orientation='horizontal')
        self.textbox_scrollbarx.grid(row=1, columnspan=2, sticky="we")

        # connect textbox scroll event to CTk scrollbar
        self.textbox.configure(yscrollcommand=self.textbox_scrollbar.set,xscrollcommand=self.textbox_scrollbarx.set)

        logo_width = 200
        self.logo_image = self.load_image("/graphics/logo.png", logo_width,int(logo_width/4.6))
        self.logo_label = customtkinter.CTkLabel(self,image=self.logo_image,width=logo_width, height=int(logo_width/4.6))
        self.logo_label.place(x=self.frame_left.winfo_rootx()+20, y=20, anchor='nw')

    #==============================================================================================================

    #==============================================================================================================
    def bootstrap_switch_samples(self):
        if self.bootstrap_switch.get():
            self.bootstrap_entry.configure(state='normal')
            self.bootstrap_entry.insert(0,'500')
        else: 
            self.bootstrap_entry.delete(0, len(self.bootstrap_entry.get()))
            self.bootstrap_entry.configure(state='disabled')
    #==============================================================================================================


    #==============================================================================================================
    def automatic_distances(self):
        rmin,rmax = dl.distancerange(self.data['t'])
        # Delete current entries
        self.rmin_entry.delete(0, len(self.rmin_entry.get()))
        self.rmax_entry.delete(0, len(self.rmax_entry.get()))
        # Enter new entries
        self.rmin_entry.insert(0,str(np.round(rmin,2)))
        self.rmax_entry.insert(0,str(np.round(rmax,2)))
    #==============================================================================================================

    #==============================================================================================================
    def load_image(self, path, image_size, image_size2=None):
        if image_size2 is None:
            image_size2 = image_size
        """ load rectangular image with path relative to PATH """
        return ImageTk.PhotoImage(Image.open(PATH + path).resize((image_size, image_size2)))
    #==============================================================================================================

    #==============================================================================================================
    def change_appearance_mode(self, new_appearance_mode):
        customtkinter.set_appearance_mode(new_appearance_mode)
    #==============================================================================================================

    #==============================================================================================================
    def on_closing(self, event=0):
        self.destroy()
    #==============================================================================================================

    #==============================================================================================================
    def generate_script(self):

        ex_model = [ex_model for ex_model in App.ex_models if ex_model in self.Exmodel_menu.get()][0] 
        bg_model = [bg_model for bg_model in App.bg_models if bg_model in self.Bmodel_menu.get()][0]      
        dd_model = [dd_model for dd_model in App.dd_models if dd_model in self.Pmodel_menu.get()][0]   
        if dd_model!='None':
            dd_model = 'dl.'+dd_model
        rmin = float(self.rmin_entry.get())            
        rmax = float(self.rmax_entry.get())            
        dr = float(self.dr_entry.get())   

        script = f"""
import numpy as np 
import deerlab as dl 
import matplotlib.pyplot as plt 

# File location
file = '{self.filepath}'

# Experimental parameters
deadtime = 0.1  # Acquisition deadtime, us \n"""
        for n in range(self.Ndelays):
            tau_value = float(getattr(self,f"delay{n+1}_entry").get())
            script += f"tau{n+1} = {tau_value}      # Inter-pulse delay #{n+1}, us \n"
        script += f"""
# Load the experimental data
t,Vexp = dl.deerload(file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t - t[0] + deadtime      # Account for deadtime

# Distance vector
r = np.arange({rmin},{rmax},{dr}) # nm

# Construct the model
experiment = dl.{ex_model}({','.join([f'tau{n+1}' for n in range(self.Ndelays)])}, pathways={[n+1 for n in range(self.Npathways) if bool(getattr(self,f'pathway{n+1}_switch').get())]})
Vmodel = dl.dipolarmodel(t,r,Pmodel={dd_model},Bmodel=dl.{bg_model}, experiment=experiment) \n"""
        
        fitargs = ['Vmodel','Vexp'] 
        regparam = self.regparam_menu.get().lower()
        fitargs.append(f'regparam={regparam}')
        if self.compactness_switch.get():
            scripte += f"dl.dipolarpenalty({dd_model},r,'compactness') \n"
            fitargs.append(f'penalties=compactness')
        if self.bootstrap_switch.get():
            bootstrap = int(self.bootstrap_entry.get())
            fitargs.append(f'bootstrap={bootstrap}')


        script += f""""
# Fit the model to the data
results = dl.fit({','.join(fitargs)})

# Print results summary
print(results)

# Extract fitted dipolar signal
Vfit = results.model
Vci = results.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)

# Plot the results
plt.figure(figsize=[6,7])
green = '{App.green}'

# Plot experimental data
plt.subplot(211)
plt.plot(t,Vexp,'.',color='grey',label='Data')

# Plot the fitted signal
plt.plot(t,Vfit,linewidth=3,label='Fit',color=green)
plt.fill_between(t,Vci[:,0],Vci[:,1],alpha=0.3,color=green)
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ ($\mu s$)')
plt.ylabel('$V(t)$ (arb.u.)')

# Plot the distance distribution
plt.subplot(212)
plt.plot(r,Pfit,linewidth=3,label='Fit',color=green)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=green,label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=green,label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (1/nm)')
plt.tight_layout()
plt.show()
        """
        file = filedialog.asksaveasfilename(defaultextension=".py")
        if file is None:
            return
        text_file = open(file, "wt", encoding='utf-8')
        n = text_file.write(script)
        text_file.close()
    #==============================================================================================================


    def openbrowser(self):
        ex_model = [ex_model for ex_model in App.ex_models if ex_model in self.Exmodel_menu.get()][0] 
        url = f"https://jeschkelab.github.io/DeerLab/_autosummary/deerlab.{ex_model}.html"
        open_new(url)

if __name__ == "__main__":
    app = App()
    app.iconbitmap(PATH + "\graphics\\favicon.ico")
    app.mainloop()