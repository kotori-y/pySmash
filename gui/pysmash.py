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


import os
import threading
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
from tkinter import Tk, Label, Entry, Button, Radiobutton
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

import pandas as pd
from getRes import getFingerprintRes


class SmashGui(Tk):
    
    def __init__(self):
        """
        Init
        
        """
        Tk.__init__(self)
        # self.pack()
        self.geometry('600x380+500+200')
        self.resizable(0,0)
        self.title('Smash molecule based on fingerprint')
        
        self.bg = '#fdd8aa'
        self.fg = '#654644'
        self.btg = '#fdafaa'
        
        self.filename = ''
        self.lblFont = ('Times New Roman', 14)
        self.createWidgets()
        try:
            self.iconbitmap('icon.ico') 
        except:
            self.iconbitmap('gui\icon.ico')  #for vscode
        
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
    
    def main_thread(self, func):
        t = threading.Thread(target=func, args=()) 
        t.setDaemon(True)
        t.start()
    
    def main(self):
        def add(words):
            textPad['state'] = 'normal'
            textPad.insert(tk.END, words)
            textPad['state'] = 'disable'
        
        self.process = tk.Toplevel(self)
        self.process.geometry('400x300+500+200')
        # self.process.resizable(0,0)
        self.process.title('Running...')
        lblnow = Label(self.process, text='Now Processing',
                       font=self.lblFont)
        lblnow.place(x=120, y=40)
        
        textPad = ScrolledText(self.process, width=40, height=15)
        textPad.place(x=55, y=80)
        textPad['state'] = 'disable'
        
        add('Load file... ')
        data = self.readFile()
        
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
        
        subMatrix, subPvalue, labels = getFingerprintRes(textPad, data, **kwgrs)

    
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
                self.cmbSmiles["values"]=self.cols
                self.cmbLabel["values"]=self.cols
                
                self.cmbSmiles['state'], self.cmbLabel['state'], self.cmbAim['state'],\
                self.cmbFP['state'] = ['readonly']*4
            else:
                disable()
            self.txtFile['state'] = 'readonly'
        
        def chooseAimLabel(*args):
            self.cmbAim.set('')
            data = self.readFile(usecols=[self.cmbLabel.get()])
            labels = list(set(data.iloc[:,0]))
            self.cmbAim['values']=labels
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
                        txtnjobs['state'], btnNext['state'] = ['disable']*11
                        
            self.cmbSmiles['state'], self.cmbLabel['state'], self.cmbAim['state'],\
                self.cmbFP['state'] = ['disable']*4            
    
        def changestate(*args):
            txtminNum['state'], txtminRatio['state'], txtPvalue['state'],\
                txtnjobs['state'], btnNext['state'] = ['normal']*5
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
        
        self.cmbSmiles = ttk.Combobox(self,width=12)
        self.cmbSmiles.place(x=85, y=100)
        
        
        lbllabel = Label(self, text='Label', 
                          font=('Times New Roman', 13),
                          bg=self.bg)
        lbllabel.place(x=210, y=100)
        
        self.cmbLabel = ttk.Combobox(self,width=12)
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
        self.cmbFP['state']="readonly"
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
        cmbSparse['state']="readonly"
        cmbSparse.bind("<<ComboboxSelected>>", ignorenBits)
        
        
        lblnBits = Label(self, text='nBits', bg=self.bg,
                          font=('Times New Roman', 13))
        lblnBits.place(x=180, y=205)
        self.nBits = tk.IntVar(value=1)
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
        
        
        ####################### Next Module#######################
        btnNext= Button(self, text='Next', 
                        font=('Times New Roman', 16),
                        bg=self.btg, width=30, height=1,
                        command=lambda: self.main_thread(self.main))
        btnNext.place(x=100, y=320)
        
        disable()
        
        
        
        
        
if '__main__' == __name__:
    gui = SmashGui()
    gui.mainloop()


