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
from getRes import getFingerprintRes, predict
from smash import ShowResult


class SmashGui(Tk):

    def __init__(self):
        """
        Init

        """
        Tk.__init__(self)
        # self.pack()
        self.geometry('600x400+500+200')
        self.resizable(0, 0)
        self.title('Smash molecule based on fingerprint')

        self.bg = '#abbfc5'
        self.fg = '#b70131'
        self.btg = '#fdafaa'

        self.filename = ''
        self.lblFont = ('Times New Roman', 14)

        self.creatTab()
        self.createWidgets()
        self.createPredictWidgets()

        self.thread_run = None
        self.thread_run_stop = Event()
        try:
            self.iconbitmap(r"icon.ico")
        except:
            self.iconbitmap(r"gui\smash\icon.ico")

    def readFile(self, file, **kwgrs):

        extendName = os.path.splitext(file)[1]
        if extendName == '.csv':
            data = pd.read_csv(file, **kwgrs)
        elif extendName == '.txt':
            data = pd.read_csv(file, sep='\t', **kwgrs)
        else:
            data = pd.read_excel(file, **kwgrs)

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

    def downloadRes(self, data, datype, preview=True, **kwgrs):
        if preview:
            self.previewRes(data, datype)
        if datype == 'df':
            savefile = asksaveasfilename(filetypes=(("CSV file", "*.csv*"), ))
            if savefile:
                try:
                    data.to_csv(savefile, **kwgrs)
                except PermissionError:
                    messagebox.showerror(
                        title='Error!', message="Permission Denied!!!")
            else:
                pass
        elif datype == 'HTML':
            savefile = asksaveasfilename(
                filetypes=(("Html file", "*.html*"), ))
            if savefile:
                try:
                    self.model.savePvalue(savefile)
                except PermissionError:
                    messagebox.showerror(
                        title='Error!', message="Permission Denied!!!")
            else:
                pass
        else:
            savefile = asksaveasfilename(
                filetypes=(("pkl file", "*.pkl*"), ))
            if savefile:
                try:
                    self.model.saveModel(savefile)
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
            self.previewPad.insert(tk.END, self.kwgrs)
            self.previewPad['state'] = 'disabled'

    def preview(self):

        self.process.destroy()
        self.view = tk.Toplevel(self)
        self.view.geometry('800x650+500+300')

        lblsubMatrix = Label(self.view, text='subMatrix',
                             fg=self.fg, font=self.lblFont)
        lblsubMatrix.place(x=100, y=80)

        lblsubPvalue = Label(self.view, text='subPvalue',
                             fg=self.fg, font=self.lblFont)
        lblsubPvalue.place(x=350, y=80)

        lblsubHTML = Label(self.view, text='Summray',
                           fg=self.fg, font=self.lblFont)
        lblsubHTML.place(x=600, y=80)

        btnPreviewMatrix = Button(self.view, text='Preview',
                                  bg=self.btg,
                                  font=('Times New Roman', 12),
                                  width=8,
                                  command=lambda: self.previewRes(
                                      data=self.subMatrix.head(), datype='df', display=5)
                                  )
        btnPreviewMatrix.place(x=60, y=120)

        btnDownloadMatrix = Button(self.view, text='Download',
                                   bg=self.btg,
                                   font=('Times New Roman', 12),
                                   width=8,
                                   command=lambda: self.downloadRes(
                                       data=self.subMatrix, datype='df', index=False)
                                   )
        btnDownloadMatrix.place(x=150, y=120)

        btnPreviewPvalue = Button(self.view, text='Preview',
                                  bg=self.btg,
                                  font=('Times New Roman', 12),
                                  width=8,
                                  command=lambda: self.previewRes(
                                      data=self.subPvalue.head(), datype='df', display=5),
                                  )
        btnPreviewPvalue.place(x=310, y=120)

        btnDownloadPvalue = Button(self.view, text='Download',
                                   bg=self.btg,
                                   font=('Times New Roman', 12),
                                   width=8,
                                   command=lambda: self.downloadRes(
                                       data=None, datype='HTML'))
        btnDownloadPvalue.place(x=400, y=120)

        btnPreviewModel = Button(self.view, text='Preview',
                                 bg=self.btg,
                                 font=('Times New Roman', 12),
                                 width=8,
                                 command=lambda: self.previewRes(
                                     None, 'Model', 50)
                                 )
        btnPreviewModel.place(x=560, y=120)

        btnDownloadModel = Button(self.view, text='Download',
                                  bg=self.btg,
                                  font=('Times New Roman', 12),
                                  width=8,
                                  command=lambda: self.downloadRes(
                                      data=None, datype='Model', escape=False))
        btnDownloadModel.place(x=650, y=120)

        self.previewPad = Text(self.view, width=105, height=35,
                               wrap="none", borderwidth=0,
                               )
        self.previewPad.place(x=20, y=160)

        vscroll = Scrollbar(self.view, orient=tk.VERTICAL,
                            command=self.previewPad.yview)
        self.previewPad['yscroll'] = vscroll.set
        vscroll.pack(side=tk.RIGHT, fill=tk.Y)

        hscroll = Scrollbar(self.view, orient=tk.HORIZONTAL,
                            command=self.previewPad.xview)
        self.previewPad['xscroll'] = hscroll.set
        hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.previewPad['state'] = 'disabled'

    def main(self):

        self.kwgrs = {'smiles_field': self.cmbSmiles.get(),
                      'label_field': self.cmbLabel.get(),
                      'fingerprint': self.cmbFP.get(),
                      'radius': self.Radius.get(),
                      'minRadius': self.minRadius.get(),
                      'minPath': self.minPath.get(),
                      'maxPath': self.maxPath.get(),
                      'folded': False,
                      'minRatio': self.minRatio.get(),
                      'minNum': self.minNum.get(),
                      'aimLabel': self.cmbAim.get(),
                      'n_jobs': self.n_jobs.get(),
                      'Bonferroni': self.Bonferroni.get(),
                      'minAccuracy': self.minAcc.get(),
                      'pValue': self.pValue.get()}

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
        data = self.readFile(self.filename)

        self.model, self.subMatrix, self.subPvalue = getFingerprintRes(
            textPad, data, **self.kwgrs)
        time.sleep(1)
        add('\nFinished!')

        btnNext['state'] = 'normal'

    def main_predict(self):

        self.detailPad['state'] = 'normal'
        self.detailPad.insert(tk.END, 'Waiting...\n')
        self.detailPad['state'] = 'disable'

        data = self.readFile(self.predFileName)
        smis = data[self.cmbPredSmiles.get()].values
        y_pred, self.predMatrix = predict(self.modelFileName, smis)
        self.predMatrix['PredLabel'] = y_pred
        self.btnSaveMatrix['state'] = 'normal'

        self.detailPad['state'] = 'normal'
        self.detailPad.insert(tk.END, 'Finished!!!\n')
        self.detailPad['state'] = 'disable'
        
        # print(self.predMatrix)

    def creatTab(self):
        tab_main = ttk.Notebook(self)
        tab_main.place(relx=0.01, rely=0.01, relwidth=0.98, relheight=0.98)

        self.fitTab = Frame(tab_main)
        tab_main.add(self.fitTab, text='Calculate')
        self.predTab = Frame(tab_main)
        tab_main.add(self.predTab, text='Predict')

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
                data = self.readFile(self.filename, nrows=0)
                self.cols = list(data.columns)
                self.cmbSmiles["values"] = self.cols
                self.cmbLabel["values"] = self.cols

                self.cmbSmiles['state'] = 'readonly'
                # self.cmbLabel['state'], self.cmbAim['state'],
                #     self.cmbFP['state'] = ['readonly']*4
            else:
                disable()
            self.txtFile['state'] = 'readonly'

        def _changesmiles(*args):
            self.cmbLabel['state'] = 'readonly'

        def chooseAimLabel(*args):
            self.cmbAim['state'] = 'readonly'
            self.cmbFP['state'] = 'readonly'
            self.cmbAim.set('')
            data = self.readFile(self.filename, usecols=[self.cmbLabel.get()])
            labels = list(set(data.iloc[:, 0]))
            self.cmbAim['values'] = labels
            self.cmbAim.current(1)

        # def ignorenBits(*args):
        #     if self.Folded.get():
        #         txtnBits['state'] = 'normal'
        #     else:
        #         txtnBits['state'] = 'disable'

        def disable(*args):
            txtminRadius['state'], txtRadius['state'], txtminPath['state'],\
                txtmaxPath['state'], txtminNum['state'], txtminRatio['state'],\
                txtPvalue['state'], txtnjobs['state'], btnRun['state'],\
                txtAcc['state'], cmbBon['state'] = ['disable']*11

            self.cmbSmiles['state'], self.cmbLabel['state'], self.cmbAim['state'],\
                self.cmbFP['state'] = ['disable']*4

        def changestate(*args):
            txtminNum['state'], txtminRatio['state'], txtPvalue['state'],\
                txtnjobs['state'], btnRun['state'], txtAcc['state'] = [
                    'normal']*6
            cmbBon['state'] = 'readonly'

            if self.cmbFP.get() == 'Circular':
                # cmbFolded['state'] = 'readonly'
                # ignorenBits()
                txtminRadius['state'], txtRadius['state'] = ['normal']*2
                txtminPath['state'], txtmaxPath['state'] = ['disable']*2

            elif self.cmbFP.get() == 'Path':
                # cmbFolded['state'] = 'normal'
                # ignorenBits()
                txtminRadius['state'], txtRadius['state'] = ['disable']*2
                txtminPath['state'], txtmaxPath['state'] = ['normal']*2

            elif self.cmbFP.get() == 'Function Group':
                txtminRadius['state'], txtRadius['state'],\
                    txtminPath['state'], txtmaxPath['state'] = ['disable']*4
        # global image
        # image = tk.PhotoImage(file='logo.gif')
        # imgLabel = Label(self, image=image).place(x=170,y=20)

        ###################### Select File Module #######################
        color = '#ffab66'
        bbg = Label(self.fitTab, bg=color,
                    width=500, height=4)
        bbg.place(x=0, y=0)

        lblFile = Label(self.fitTab, text='>>> Select the file',
                        font=self.lblFont, bg=color,
                        fg='#b70131')
        lblFile.place(x=5, y=10)

        self.txtFile = Entry(self.fitTab, width=60)
        self.txtFile.place(x=7, y=35)
        self.txtFile['state'] = 'readonly'

        btnGetFile = Button(self.fitTab, text='Browse...',
                            command=getFileName,
                            bg='#66baff',
                            width=18)
        btnGetFile.place(x=440, y=30)
        ####################### Select File Module #######################

        ####################### Select Aim Field Module #######################
        color = '#ffb97f'
        bbg = Label(self.fitTab, bg=color,
                    width=500, height=4)
        bbg.place(x=0, y=74)

        lblField = Label(self.fitTab, text='>>> Select related field',
                         font=self.lblFont, bg=color, fg=self.fg)
        lblField.place(x=0, y=74)

        lblSmiles = Label(self.fitTab, text='SMILES',
                          font=('Times New Roman', 12),
                          bg=color)
        lblSmiles.place(x=20, y=105)

        self.cmbSmiles = ttk.Combobox(self.fitTab, width=12)
        self.cmbSmiles.place(x=85, y=105)

        self.cmbSmiles.bind("<<ComboboxSelected>>", _changesmiles)

        lbllabel = Label(self.fitTab, text='Label',
                         font=('Times New Roman', 13),
                         bg=color)
        lbllabel.place(x=210, y=105)

        self.cmbLabel = ttk.Combobox(self.fitTab, width=12)
        self.cmbLabel.place(x=260, y=105)
        self.cmbLabel.bind("<<ComboboxSelected>>", chooseAimLabel)

        lbllabel = Label(self.fitTab, text='Aim Label',
                         font=('Times New Roman', 13),
                         bg=color)
        lbllabel.place(x=385, y=105)

        self.cmbAim = ttk.Combobox(self.fitTab, width=12)
        self.cmbAim.place(x=468, y=105)
        ####################### Select Aim Field Module #######################

        ####################### Select Fragment Type #######################
        color = '#ffc799'
        bbg = Label(self.fitTab, bg=color,
                    width=45, height=10)
        bbg.place(x=0, y=140)

        lblFPM = Label(self.fitTab, text=">>> Adjust fragment parameter",
                       font=self.lblFont, bg=color, fg=self.fg)
        lblFPM.place(x=0, y=140)

        lblFP = Label(self.fitTab, text='Fragment Type',
                      font=('Times New Roman', 12),
                      bg=color)
        lblFP.place(x=15, y=180)

        self.cmbFP = ttk.Combobox(self.fitTab, width=14)
        self.cmbFP['values'] = ['Circular', 'Path', 'Function Group']
        self.cmbFP.place(x=120, y=180)
        self.cmbFP['state'] = "readonly"
        self.cmbFP.bind("<<ComboboxSelected>>", changestate)
        ####################### Select Fragment Type #######################

        ####################### Adjust Figerprint Param  Module#######################
        # lblFolded = Label(self, text='Fold', bg=self.bg,
        #                   font=('Times New Roman', 13))
        # lblFolded.place(x=40, y=205)
        # self.Folded = tk.BooleanVar()
        # cmbFolded = ttk.Combobox(self, width=5, textvariable=self.Folded)
        # cmbFolded['values'] = [True, False]
        # cmbFolded.current(1)
        # cmbFolded.place(x=103, y=205)
        # cmbFolded['state'] = "readonly"
        # cmbFolded.bind("<<ComboboxSelected>>", ignorenBits)

        # lblnBits = Label(self, text='nBits', bg=self.bg,
        #                  font=('Times New Roman', 13))
        # lblnBits.place(x=180, y=205)
        # self.nBits = tk.IntVar(value=1024)
        # txtnBits = Entry(self, width=6, textvariable=self.nBits)
        # txtnBits.place(x=243, y=205)

        lblminRadius = Label(self.fitTab, text='minRadius', bg=color,
                             font=('Times New Roman', 13))
        lblminRadius.place(x=15, y=220)
        self.minRadius = tk.IntVar(value=1)
        txtminRadius = Entry(self.fitTab, width=5, textvariable=self.minRadius)
        txtminRadius.place(x=95, y=220)

        lblRadius = Label(self.fitTab, text='maxRadius', bg=color,
                          font=('Times New Roman', 13))
        lblRadius.place(x=155, y=220)
        self.Radius = tk.IntVar(value=2)
        txtRadius = Entry(self.fitTab, width=5, textvariable=self.Radius)
        txtRadius.place(x=235, y=220)

        lblminPath = Label(self.fitTab, text='minPath', bg=color,
                           font=('Times New Roman', 13))
        lblminPath.place(x=15, y=275)
        self.minPath = tk.IntVar(value=1)
        txtminPath = Entry(self.fitTab, width=5, textvariable=self.minPath)
        txtminPath.place(x=95, y=275)

        lblmaxPath = Label(self.fitTab, text='maxPath', bg=color,
                           font=('Times New Roman', 13))
        lblmaxPath.place(x=155, y=275)
        self.maxPath = tk.IntVar(value=7)
        txtmaxPath = Entry(self.fitTab, width=5, textvariable=self.maxPath)
        txtmaxPath.place(x=235, y=275)
        ####################### Adjust Figerprint Param Module#######################

        ####################### Adjust Running Param Module#######################
        color = '#ffd5b2'
        bbg = Label(self.fitTab, bg=color,
                    width=45, height=10)
        bbg.place(x=310, y=140)
        lblRP = Label(self.fitTab, text='>>> Adjust running parameter',
                      bg=color, fg=self.fg, font=self.lblFont)
        lblRP.place(x=310, y=140)

        lblminNum = Label(self.fitTab, text='minNum', bg=color,
                          font=('Times New Roman', 13))
        lblminNum.place(x=320, y=180)
        self.minNum = tk.IntVar(value=5)
        txtminNum = Entry(self.fitTab, width=7, textvariable=self.minNum)
        txtminNum.place(x=390, y=180)

        lblRatio = Label(self.fitTab, text='minRatio', bg=color,
                         font=('Times New Roman', 13))
        lblRatio.place(x=450, y=180)
        self.minRatio = tk.DoubleVar(value=0.4)
        txtminRatio = Entry(self.fitTab, width=7, textvariable=self.minRatio)
        txtminRatio.place(x=520, y=180)

        lblPvalue = Label(self.fitTab, text='p-value', bg=color,
                          font=('Times New Roman', 13))
        lblPvalue.place(x=320, y=230)
        self.pValue = tk.DoubleVar(value=0.05)
        txtPvalue = Entry(self.fitTab, width=7, textvariable=self.pValue)
        txtPvalue.place(x=390, y=230)

        lblAcc = Label(self.fitTab, text='minAcc', bg=color,
                       font=('Times New Roman', 13))
        lblAcc.place(x=450, y=230)
        self.minAcc = tk.DoubleVar(value=0.70)
        txtAcc = Entry(self.fitTab, width=7, textvariable=self.minAcc)
        txtAcc.place(x=520, y=230)

        lblnjobs = Label(self.fitTab, text='n_jobs', bg=color,
                         font=('Times New Roman', 13))
        lblnjobs.place(x=320, y=280)
        self.n_jobs = tk.IntVar(value=1)
        txtnjobs = Entry(self.fitTab, width=7, textvariable=self.n_jobs)
        txtnjobs.place(x=390, y=280)

        lblBon = Label(self.fitTab, text='Bonferroni',
                       font=('Times New Roman', 12),
                       bg=color)
        lblBon.place(x=450, y=280)

        self.Bonferroni = tk.BooleanVar()
        cmbBon = ttk.Combobox(self.fitTab, width=4,
                              textvariable=self.Bonferroni)
        cmbBon['values'] = [False, True]
        cmbBon.current(0)
        cmbBon.place(x=520, y=280)

        ####################### Adjust Running Param Module#######################

        ####################### Run Module#######################
        color = '#fff1e5'
        bbg = Label(self.fitTab, bg=color,
                    width=100, height=10)
        bbg.place(x=0, y=310)

        btnRun = Button(self.fitTab, text='Calculate',
                        font=('Times New Roman', 16),
                        bg='#e5f3ff', width=10, height=1,
                        command=lambda: self.main_thread(self.main),
                        # command=self.preview
                        )
        btnRun.place(x=210, y=320)

        disable()

    def createPredictWidgets(self):
        #####################################################################
        # Prediction
        # Prediction
        # Prediction
        # Prediction
        #####################################################################
        def getFileName():
            self.txtPredFile['state'] = 'normal'
            self.txtPredFile.delete(0, tk.END)
            self.predFileName = askopenfilename(
                filetypes=(("csv file", "*.csv*"),
                           ("Excel file", "*.xlsx*;*.xls*"),
                           ("Text file", "*.txt*")))
            if self.predFileName:
                self.txtPredFile.insert(tk.END, self.predFileName)
                data = self.readFile(self.predFileName, nrows=0)
                self.predCols = list(data.columns)
                self.cmbPredSmiles["values"] = self.predCols

                self.cmbPredSmiles['state'] = 'readonly'
            # else:
            #     disable()
            self.txtPredFile['state'] = 'readonly'

        def getModelFileName():
            self.txtModelFile['state'] = 'normal'
            self.txtModelFile.delete(0, tk.END)
            self.modelFileName = askopenfilename(
                filetypes=(("pbz2 file", "*.pbz2*"),))
            if self.modelFileName:
                self.txtModelFile.insert(tk.END, self.modelFileName)
                # data = self.readFile(self.PvFileName, nrows=0)
                # self.pvCols = list(data.columns)
                # self.cmbPvalue["values"] = self.pvCols

                # self.cmbPvalue['state'] = 'readonly'
            # else:
            #     disable()
            self.txtModelFile['state'] = 'readonly'

        color = '#ffab66'
        bbg = Label(self.predTab, bg=color,
                    width=500, height=4)
        bbg.place(x=0, y=0)

        lblPredFile = Label(self.predTab, text='>>> Select the file and SMILES field',
                            font=self.lblFont, bg=color,
                            fg='#b70131')
        lblPredFile.place(x=0, y=10)

        self.txtPredFile = Entry(self.predTab, width=50)
        self.txtPredFile.place(x=7, y=35)
        self.txtPredFile['state'] = 'readonly'

        btnGetPredFile = Button(self.predTab, text='Browse...',
                                command=getFileName,
                                bg='#66baff',
                                width=7)
        btnGetPredFile.place(x=365, y=30)

        self.cmbPredSmiles = ttk.Combobox(self.predTab, width=12)
        self.cmbPredSmiles.place(x=450, y=30)
        self.cmbPredSmiles['state'] = 'disable'

        color = '#ffb97f'
        bbg = Label(self.predTab, bg=color,
                    width=500, height=4)
        bbg.place(x=0, y=74)

        lblModelFile = Label(self.predTab, text='>>> Select the model file',
                             font=self.lblFont, bg=color,
                             fg='#b70131')
        lblModelFile.place(x=0, y=74)

        self.txtModelFile = Entry(self.predTab, width=70)
        self.txtModelFile.place(x=7, y=110)
        self.txtModelFile['state'] = 'readonly'

        btnGetModelFile = Button(self.predTab, text='Browse...',
                                 command=getModelFileName,
                                 bg='#66baff',
                                 width=7)
        btnGetModelFile.place(x=505, y=105)

        # self.cmbPvalue = ttk.Combobox(self.predTab, width=12)'
        # self.cmbPvalue.place(x=450, y=110)
        # self.cmbPvalue['state'] = 'disable'

        color = '#ffd5b2'
        bbg = Label(self.predTab, bg=color,
                    width=500, height=7)
        bbg.place(x=0, y=147)

        btnPredict = Button(self.predTab, text='Predict',
                            command=lambda: self.main_thread(
                                self.main_predict),
                            bg='#66baff',
                            width=7)
        btnPredict.place(x=250, y=165)

        # self.btnSaveLabel = Button(self.predTab, text='Save Predict Label',
        #                            command=lambda: self.main_thread(
        #                                self.main_predict),
        #                            bg='#66baff',
        #                            width=20)
        # self.btnSaveLabel.place(x=140, y=220)
        # self.btnSaveLabel['state'] = 'disable'

        self.btnSaveMatrix = Button(self.predTab, text='Save Predict Result',
                                    command=lambda: self.downloadRes(
                                        data=self.predMatrix, datype='df', preview=False, index=False),
                                    bg='#66baff',
                                    width=20)
        self.btnSaveMatrix.place(x=210, y=220)
        self.btnSaveMatrix['state'] = 'disable'

        self.detailPad = Text(self.predTab, width=30, height=5,
                              wrap="none", borderwidth=0,
                              )
        self.detailPad.place(x=170, y=280)


if '__main__' == __name__:

    mp.freeze_support()
    gui = SmashGui()
    gui.mainloop()
