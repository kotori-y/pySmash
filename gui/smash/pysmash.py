# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 17:15:00 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


import multiprocessing as mp
import os
import time
import socket
from threading import Thread, Event
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showinfo
from tkinter.scrolledtext import ScrolledText
from tkinter import Tk, Label, Entry, Button, Radiobutton, Scrollbar, Text, Frame
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter import messagebox

import pandas as pd
from getRes import getFingerprintRes
from smash import ShowResult


class SmashGui(Tk):

    def __init__(self):
        """
        Init

        """
        Tk.__init__(self)
        # self.pack()
        self.geometry('600x380+500+200')
        self.resizable(0, 0)
        self.title('Smash molecule based on fingerprint')

        self.bg = '#fdd8aa'
        self.fg = '#654644'
        self.btg = '#fdafaa'

        self.filename = ''
        self.lblFont = ('Times New Roman', 14)
        self.createWidgets()

        self.thread_run = None
        self.thread_run_stop = Event()
        try:
            self.iconbitmap(r"icon.ico")
        except:
            self.iconbitmap(r"gui\smash\icon.ico")

    def readFile(self, **kwgrs):

        extendName = os.path.splitext(self.filename)[1]
        if extendName == '.csv':
            data = pd.read_csv(self.filename, **kwgrs)
        elif extendName == '.txt':
            data = pd.read_csv(self.filename, sep='\t', **kwgrs)
        else:
            data = pd.read_excel(self.filename, **kwgrs)

        return data
        # print(self.data)

    def main_thread(self, func, args=()):
        self.thread_run = Thread(target=func, args=args)
        self.thread_run.setDaemon(True)
        self.thread_run.start()

    def stop_run(self):
        self.thread_run_stop.set()
        self.thread_run.join()
        self.thread_run = None

    def downloadRes(self, data, datype, **kwgrs):
        self.previewRes(data, datype)
        if datype == 'df':
            savefile = asksaveasfilename(defaultextension=".csv")
            if savefile:
                try:
                    data.to_csv(savefile, **kwgrs)
                except PermissionError:
                    messagebox.showerror(
                        title='Error!', message="Permission Denied!!!")
            else:
                pass
        else:
            savefile = asksaveasfilename(defaultextension=".html")
            if savefile:
                try:
                    self.html.to_html(savefile, **kwgrs)
                except PermissionError:
                    messagebox.showerror(
                        title='Error!', message="Permission Denied!!!")
            else:
                pass

    def previewRes(self, data, datype=None, display=5):
        self.previewPad['state'] = 'normal'
        self.previewPad.delete(1.0, tk.END)
        if datype == 'df':
            data = data.head(display)
            self.previewPad.insert(tk.END, data.to_string(
                max_cols=10, justify='center'))
            self.previewPad['state'] = 'disabled'
        else:
            try:
                aimLabel = float(self.cmbAim.get())
            except:
                pass
            self.html = ShowResult(self.subMatrix, self.subPvalue,
                                   smarts_field=None, aimLabel=aimLabel,
                                   topx=display)

            self.previewPad.insert(
                tk.END, self.html.iloc[:, :-1].to_string(max_cols=10, justify='center'))
            self.previewPad['state'] = 'disabled'

    def preview(self):

        view = tk.Toplevel(self)
        view.geometry('800x650+500+300')

        lblsubMatrix = Label(view, text='subMatrix',
                             fg=self.fg, font=self.lblFont)
        lblsubMatrix.place(x=100, y=80)

        lblsubPvalue = Label(view, text='subPvalue',
                             fg=self.fg, font=self.lblFont)
        lblsubPvalue.place(x=350, y=80)

        lblsubHTML = Label(view, text='Summray',
                           fg=self.fg, font=self.lblFont)
        lblsubHTML.place(x=600, y=80)

        btnPreviewMatrix = Button(view, text='Preview',
                                  bg=self.btg,
                                  font=('Times New Roman', 12),
                                  width=8,
                                  command=lambda: self.previewRes(
                                      data=self.subMatrix.head(), datype='df', display=5)
                                  )
        btnPreviewMatrix.place(x=60, y=120)

        btnDownloadMatrix = Button(view, text='Download',
                                   bg=self.btg,
                                   font=('Times New Roman', 12),
                                   width=8,
                                   command=lambda: self.downloadRes(
                                       data=self.subMatrix, datype='df', index=False)
                                   )
        btnDownloadMatrix.place(x=150, y=120)

        btnPreviewPvalue = Button(view, text='Preview',
                                  bg=self.btg,
                                  font=('Times New Roman', 12),
                                  width=8,
                                  command=lambda: self.previewRes(
                                      data=self.subPvalue.head(), datype='df', display=5),
                                  )
        btnPreviewPvalue.place(x=310, y=120)

        btnDownloadPvalue = Button(view, text='Download',
                                   bg=self.btg,
                                   font=('Times New Roman', 12),
                                   width=8,
                                   command=lambda: self.downloadRes(
                                       data=self.subPvalue, datype='df', index_label='SMARTS'))
        btnDownloadPvalue.place(x=400, y=120)

        btnPreviewHTML = Button(view, text='Preview',
                                bg=self.btg,
                                font=('Times New Roman', 12),
                                width=8,
                                command=lambda: self.previewRes(
                                    None, 'HTML', 50)
                                )
        btnPreviewHTML.place(x=560, y=120)

        btnDownloadHTML = Button(view, text='Download',
                                 bg=self.btg,
                                 font=('Times New Roman', 12),
                                 width=8,
                                 command=lambda: self.downloadRes(
                                     data=None, datype='HTML', escape=False))
        btnDownloadHTML.place(x=650, y=120)

        self.previewPad = Text(view, width=105, height=35,
                               wrap="none", borderwidth=0,
                               )
        self.previewPad.place(x=20, y=160)

        vscroll = Scrollbar(view, orient=tk.VERTICAL,
                            command=self.previewPad.yview)
        self.previewPad['yscroll'] = vscroll.set
        vscroll.pack(side=tk.RIGHT, fill=tk.Y)

        hscroll = Scrollbar(view, orient=tk.HORIZONTAL,
                            command=self.previewPad.xview)
        self.previewPad['xscroll'] = hscroll.set
        hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.previewPad['state'] = 'disabled'

    def main(self):

        kwgrs = {'smiles_field': self.cmbSmiles.get(),
                    'label_field': self.cmbLabel.get(),
                    'fingerprint': self.cmbFP.get(),
                    'radius': self.Radius.get(),
                    'minRadius': self.minRadius.get(),
                    'minPath': self.minPath.get(),
                    'maxPath': self.maxPath.get(),
                    'nBits': self.nBits.get(),
                    'sparse': self.Sparse.get(),
                    'minRatio': self.minRatio.get(),
                    'minNum': self.minNum.get(),
                    'aimLabel': self.cmbAim.get(),
                    'n_jobs': self.n_jobs.get()}
        


        def add(words):
            textPad['state'] = 'normal'
            textPad.insert(tk.END, words)
            textPad['state'] = 'disable'

        self.process = tk.Toplevel(self)
        self.process.geometry('400x300+500+200')
        # self.process.resizable(0,0)
        self.process.title('Running...')
        lblnow = Label(self.process, text='Processing',
                       font=self.lblFont)
        lblnow.place(x=135, y=40)

        textPad = ScrolledText(self.process, width=48, height=13)
        textPad.place(x=43, y=85)
        textPad['state'] = 'disable'

        btnNext = Button(self.process, text='Next', command=self.preview)
        btnNext.place(x=320, y=265, width=50, height=25)
        btnNext['state'] = 'disable'

        btnCancel = Button(self.process, text='Cancel',
                           command=lambda: self.process.destroy())
        btnCancel.place(x=260, y=265, width=50, height=25)

        add('Load file... ')
        data = self.readFile()

        self.subMatrix, self.subPvalue = getFingerprintRes(textPad,
                                                        data, **kwgrs)
        time.sleep(1)
        add('\nFinished!')

        btnNext['state'] = 'normal'

    def createWidgets(self):
        def getFileName():
            self.txtFile['state'] = 'normal'
            self.txtFile.delete(0, tk.END)
            self.filename = askopenfilename(
                filetypes=(("csv file", "*.csv*"),
                           ("Excel file", "*.xlsx*;*.xls*"),
                           ("Text file", "*.txt*")))
            if self.filename:
                self.txtFile.insert(tk.END, self.filename)
                data = self.readFile(nrows=0)
                self.cols = list(data.columns)
                self.cmbSmiles["values"] = self.cols
                self.cmbLabel["values"] = self.cols

                self.cmbSmiles['state'] = 'readonly'
                # self.cmbLabel['state'], self.cmbAim['state'],\
                #     self.cmbFP['state'] = ['readonly']*4
            else:
                disable()
            self.txtFile['state'] = 'readonly'

        def _changesmiles(*args):
                    self.cmbLabel['state']='readonly'

        def chooseAimLabel(*args):
            self.cmbAim['state']='readonly'
            self.cmbFP['state']='readonly'
            self.cmbAim.set('')
            data = self.readFile(usecols=[self.cmbLabel.get()])
            labels = list(set(data.iloc[:, 0]))
            self.cmbAim['values'] = labels
            self.cmbAim.current(1)

        def ignorenBits(*args):
            if self.Sparse.get():
                txtnBits['state'] = 'disable'
            else:
                txtnBits['state'] = 'normal'

        def disable(*args):
            cmbSparse['state'], txtnBits['state'], txtminRadius['state'],\
                txtRadius['state'], txtminPath['state'], txtmaxPath['state'],\
                txtminNum['state'], txtminRatio['state'], txtPvalue['state'],\
                txtnjobs['state'], btnRun['state'] = ['disable']*11

            self.cmbSmiles['state'], self.cmbLabel['state'], self.cmbAim['state'],\
                self.cmbFP['state'] = ['disable']*4

        def changestate(*args):
            txtminNum['state'], txtminRatio['state'], txtPvalue['state'],\
                txtnjobs['state'], btnRun['state'] = ['normal']*5
            if self.cmbFP.get() == 'ECFP':
                cmbSparse['state'] = 'normal'
                ignorenBits()
                txtminRadius['state'], txtRadius['state'] = ['normal']*2
                txtminPath['state'], txtmaxPath['state'] = ['disable']*2

            elif self.cmbFP.get() == 'Daylight':
                cmbSparse['state'] = 'normal'
                ignorenBits()
                txtminRadius['state'], txtRadius['state'] = ['disable']*2
                txtminPath['state'], txtmaxPath['state'] = ['normal']*2

        # global image
        # image = tk.PhotoImage(file='logo.gif')
        # imgLabel = Label(self, image=image).place(x=170,y=20)

        bbg = Label(self, bg=self.bg,
                    width=500, height=300)
        bbg.pack()

        ###################### Select File Module #######################
        lblFile = Label(self, text='Select a file:',
                        font=self.lblFont, bg=self.bg,
                        fg=self.fg)
        lblFile.place(x=5, y=10)

        self.txtFile = Entry(self, width=60)
        self.txtFile.place(x=7, y=35)
        self.txtFile['state'] = 'readonly'

        btnGetFile = Button(self, text='Browse...',
                            command=getFileName,
                            bg=self.btg,
                            width=20)
        btnGetFile.place(x=440, y=30)

        ####################### Select File Module #######################

        ####################### Select Aim Field Module #######################
        lblField = Label(self, text='Select related field and aim label:',
                         font=self.lblFont, bg=self.bg, fg=self.fg)
        lblField.place(x=7, y=70)

        lblSmiles = Label(self, text='SMILES',
                          font=('Times New Roman', 12),
                          bg=self.bg)
        lblSmiles.place(x=20, y=100)

        self.cmbSmiles = ttk.Combobox(self, width=12)
        self.cmbSmiles.place(x=85, y=100)
        
        self.cmbSmiles.bind("<<ComboboxSelected>>", _changesmiles)

        lbllabel = Label(self, text='Label',
                         font=('Times New Roman', 13),
                         bg=self.bg)
        lbllabel.place(x=210, y=100)

        self.cmbLabel = ttk.Combobox(self, width=12)
        self.cmbLabel.place(x=260, y=100)
        self.cmbLabel.bind("<<ComboboxSelected>>", chooseAimLabel)

        lbllabel = Label(self, text='Aim Label',
                         font=('Times New Roman', 13),
                         bg=self.bg)
        lbllabel.place(x=385, y=100)

        self.cmbAim = ttk.Combobox(self, width=12)
        self.cmbAim.place(x=468, y=100)
        ####################### Select Aim Field Module #######################

        ####################### Select Fingerprint Module #######################
        lblFPM = Label(self, text="Select Fingerprint and Adjust Param",
                       font=self.lblFont, bg=self.bg, fg=self.fg)
        lblFPM.place(x=27, y=140)

        lblFP = Label(self, text='Fingerprint',
                      font=('Times New Roman', 13),
                      bg=self.bg)
        lblFP.place(x=50, y=170)

        self.cmbFP = ttk.Combobox(self, width=14)
        self.cmbFP['values'] = ['ECFP', 'Daylight']
        self.cmbFP.place(x=135, y=170)
        self.cmbFP['state'] = "readonly"
        self.cmbFP.bind("<<ComboboxSelected>>", changestate)
        ####################### Select Fingerprint Module #######################

        ####################### Adjust Figerprint Param  Module#######################
        lblSparse = Label(self, text='Sparse', bg=self.bg,
                          font=('Times New Roman', 13))
        lblSparse.place(x=40, y=205)
        self.Sparse = tk.BooleanVar()
        cmbSparse = ttk.Combobox(self, width=5, textvariable=self.Sparse)
        cmbSparse['values'] = [True, False]
        cmbSparse.current(0)
        cmbSparse.place(x=103, y=205)
        cmbSparse['state'] = "readonly"
        cmbSparse.bind("<<ComboboxSelected>>", ignorenBits)

        lblnBits = Label(self, text='nBits', bg=self.bg,
                         font=('Times New Roman', 13))
        lblnBits.place(x=180, y=205)
        self.nBits = tk.IntVar(value=1024)
        txtnBits = Entry(self, width=6, textvariable=self.nBits)
        txtnBits.place(x=243, y=205)

        lblminRadius = Label(self, text='minRadius', bg=self.bg,
                             font=('Times New Roman', 13))
        lblminRadius.place(x=40, y=240)
        self.minRadius = tk.IntVar(value=1)
        txtminRadius = Entry(self, width=5, textvariable=self.minRadius)
        txtminRadius.place(x=120, y=240)

        lblRadius = Label(self, text='Radius', bg=self.bg,
                          font=('Times New Roman', 13))
        lblRadius.place(x=180, y=240)
        self.Radius = tk.IntVar(value=2)
        txtRadius = Entry(self, width=5, textvariable=self.Radius)
        txtRadius.place(x=249, y=240)

        lblminPath = Label(self, text='minPath', bg=self.bg,
                           font=('Times New Roman', 13))
        lblminPath.place(x=40, y=275)
        self.minPath = tk.IntVar(value=1)
        txtminPath = Entry(self, width=5, textvariable=self.minPath)
        txtminPath.place(x=120, y=275)

        lblmaxPath = Label(self, text='maxPath', bg=self.bg,
                           font=('Times New Roman', 13))
        lblmaxPath.place(x=180, y=275)
        self.maxPath = tk.IntVar(value=7)
        txtmaxPath = Entry(self, width=5, textvariable=self.maxPath)
        txtmaxPath.place(x=250, y=275)

        ####################### Adjust Figerprint Param Module#######################

        ####################### Adjust Running Param Module#######################
        lblRP = Label(self, text='Adjust Running Param:',
                      bg=self.bg, fg=self.fg, font=self.lblFont)
        lblRP.place(x=340, y=140)

        lblminNum = Label(self, text='minNum', bg=self.bg,
                          font=('Times New Roman', 13))
        lblminNum.place(x=360, y=170)
        self.minNum = tk.IntVar(value=5)
        txtminNum = Entry(self, width=7, textvariable=self.minNum)
        txtminNum.place(x=440, y=170)

        lblRatio = Label(self, text='minRatio', bg=self.bg,
                         font=('Times New Roman', 13))
        lblRatio.place(x=360, y=205)
        self.minRatio = tk.DoubleVar(value=0.4)
        txtminRatio = Entry(self, width=7, textvariable=self.minRatio)
        txtminRatio.place(x=440, y=205)

        lblPvalue = Label(self, text='p-value', bg=self.bg,
                          font=('Times New Roman', 13))
        lblPvalue.place(x=360, y=240)
        self.Pvalue = tk.DoubleVar(value=0.05)
        txtPvalue = Entry(self, width=7, textvariable=self.Pvalue)
        txtPvalue.place(x=440, y=240)

        lblnjobs = Label(self, text='n_jobs', bg=self.bg,
                         font=('Times New Roman', 13))
        lblnjobs.place(x=360, y=275)
        self.n_jobs = tk.IntVar(value=1)
        txtnjobs = Entry(self, width=7, textvariable=self.n_jobs)
        txtnjobs.place(x=440, y=275)
        ####################### Adjust Running Param Module#######################

        ####################### Run Module#######################
        btnRun = Button(self, text='Run',
                        font=('Times New Roman', 16),
                        bg=self.btg, width=30, height=1,
                        command=lambda: self.main_thread(self.main),
                        # command=self.preview
                        )
        btnRun.place(x=100, y=320)

        disable()


if '__main__' == __name__:
    mp.freeze_support()
    gui = SmashGui()
    gui.mainloop()
