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
import tkinter as tk
from tkinter import ttk
from tkinter import Tk, Label, Entry, Button, Radiobutton
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

import pandas as pd



class SmashGui(Tk):
    
    def __init__(self):
        """
        Init
        
        """
        Tk.__init__(self)
        # self.pack()
        self.geometry('600x300+500+200')
        self.title('Smash molecule based on fingerprint')
        
        self.bg = '#fdd8aa'
        self.fg = '#654644'
        self.btg = '#fdafaa'
        
        self.filename = ''
        self.lblFont = ('Times New Roman', 14)
        self.createWidgets()
        self.iconbitmap('icon.ico')  
        
    def getFileName(self):
        
        self.filename = askopenfilename(
            filetypes=(("csv file", "*.csv*"), 
                       ("Excel file", "*.xlsx*;*.xls*"), 
                       ("Text file", "*.txt*")))
        
        self.txtFile['state'] = 'normal'
        self.txtFile.delete(0, tk.END)
        self.txtFile.insert(tk.END, self.filename)
        self.txtFile['state'] = 'readonly'
        
        data = self.readFile(nrows=0)
        self.cols = list(data.columns)
        self.cmbSmiles["values"]=self.cols
        self.cmbLabel["values"]=self.cols
        # data = self.readFile(usecols=)
    
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
           
    def createWidgets(self):
        def chooseAimLabel(*args):
            cmbAim.set('')
            data = self.readFile(usecols=[self.cmbLabel.get()])
            labels = list(set(data.iloc[:,0]))
            cmbAim['values']=labels
            cmbAim.current(1)
            
        def ignorenBits(*args):
            if self.Sparse.get():
                cmbnBits.delete(0, tk.END)
                cmbnBits.insert(tk.END, '---')
                cmbnBits['state'] = 'disable'
            else:
                cmbnBits['state'] = 'normal'
                cmbnBits.delete(0, tk.END)
                cmbnBits.insert(tk.END, 2048)
                
            
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
                            command=self.getFileName,
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
        self.cmbSmiles['state']="readonly"
        
        
        lbllabel = Label(self, text='Label', 
                          font=('Times New Roman', 13),
                          bg=self.bg)
        lbllabel.place(x=210, y=100)
        
        self.cmbLabel = ttk.Combobox(self,width=12)
        self.cmbLabel.place(x=260, y=100)
        self.cmbLabel['state']="readonly"
        self.cmbLabel.bind("<<ComboboxSelected>>", chooseAimLabel)
        
        
        lbllabel = Label(self, text='Aim Label', 
                          font=('Times New Roman', 13),
                          bg=self.bg)
        lbllabel.place(x=385, y=100)
        
        self.aimLabel = tk.IntVar()
        cmbAim = ttk.Combobox(self, width=12, textvariable=self.aimLabel)
        cmbAim.place(x=468, y=100)
        cmbAim['state']="readonly"
        ####################### Select Aim Field Module #######################
        
        
        ####################### Select Fingerprint Module #######################
        lblFP = Label(self, text="Select fingerprint and adjust param:",
                      font=self.lblFont, bg=self.bg, fg=self.fg)
        lblFP.place(x=7, y=140)
        
        lblSmiles = Label(self, text='Fingerprint', 
                          font=('Times New Roman', 13),
                          bg=self.bg)
        lblSmiles.place(x=20, y=170)
        
        self.cmbFP = ttk.Combobox(self, width=14)
        self.cmbFP['values'] = ['ECFP', 'Daylight']
        self.cmbFP.place(x=110, y=170)
        self.cmbFP['state']="readonly"
        ####################### Select Fingerprint Module #######################
        
        
        ####################### Adjust Figerprint Param  Module#######################
        lblSparse = Label(self, text='Sparse', bg=self.bg,
                          font=('Times New Roman', 13))
        lblSparse.place(x=20, y=205)
        
        self.Sparse = tk.BooleanVar()
        cmbSparse = ttk.Combobox(self, width=7, textvariable=self.Sparse)
        cmbSparse['values'] = [True, False]
        cmbSparse.current(0)
        cmbSparse.place(x=70, y=205)
        cmbSparse['state']="readonly"
        cmbSparse.bind("<<ComboboxSelected>>", ignorenBits)
        
        
        lblnBits = Label(self, text='nBits', bg=self.bg,
                          font=('Times New Roman', 13))
        lblnBits.place(x=160, y=205)
        
        self.nBits = tk.IntVar()
        cmbnBits = Entry(self, width=7, textvariable=self.nBits)
        cmbnBits.place(x=205, y=205)
        if self.Sparse.get():
            cmbnBits.delete(0, tk.END)
            cmbnBits.insert(tk.END, '---')
            cmbnBits['state']="disable"
        
        
        

        
if '__main__' == __name__:
    gui = SmashGui()
    gui.mainloop()
    # gui.readFile()