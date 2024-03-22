import pandas as pd
from rdkit import Chem
from molvs import standardize_smiles
from rdkit.Chem import SaltRemover
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
from PIL import ImageTk, Image
#import pymysql
import os
import shutil
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askopenfilenames
import warnings
#import cirpy

warnings.simplefilter(action='ignore')
initialdir=os.getcwd()

def data1():
    filename = askopenfilename(initialdir=initialdir,title = "Select excel file")
    T1E1.delete(0, END)
    T1E1.insert(0, filename)
    global a_
    a_,b_=os.path.splitext(filename)
    global file1
    file1 = pd.read_csv(filename)

def data2():
    filename = askopenfilename(initialdir=initialdir,title = "Select first input file")
    T1E2.delete(0, END)
    T1E2.insert(0, filename)
    global a_
    a_,b_=os.path.splitext(filename)
    global file2
    file2 = pd.read_csv(filename)

def data3():
    filename = askopenfilename(initialdir=initialdir,title = "Upload modified ECOTOX file")
    T2E1.delete(0, END)
    T2E1.insert(0, filename)
    global a_
    a_,b_=os.path.splitext(filename)
    global file3
    file3 = pd.read_csv(filename)

def data4():
    filename = askopenfilename(initialdir=initialdir,title = "Upload SMILES file")
    T2E2.delete(0, END)
    T2E2.insert(0, filename)
    global a_
    a_,b_=os.path.splitext(filename)
    global file4
    file4 = pd.read_csv(filename, encoding = "ISO-8859-1")

def data5():
    filename = askopenfilename(initialdir=initialdir,title = "Upload second input file")
    T2E3.delete(0, END)
    T2E3.insert(0, filename)
    global a_
    a_,b_=os.path.splitext(filename)
    global file5
    file5 = pd.read_csv(filename)


#Select columns based on first input
def process(df,inp,l1):
    m,n=[],[]
    for i in range(len(l1)):
        m.append(inp[inp[l1[i]]!='NA'][l1[i]].tolist())
        n.append(l1[i])
    for i,j in zip(m,n):
        d1=df[df[j].isin(i)]
        df=d1
    return df[inp.columns]


def castosmi(cas):
    hcas=cas[0:-3]+'-'+cas[-3:-1]+'-'+cas[-1:]
    smi=cirpy.resolve(hcas, 'smiles')
    return smi

remover = SaltRemover.SaltRemover()
def processSmi(smile):
    smile=standardize_smiles(smile)
    m=Chem.MolFromSmiles(smile)
    ms = remover.StripMol(m,dontRemoveEverything=True)
    largest_Fragment = rdMolStandardize.LargestFragmentChooser()
    largest_mol = largest_Fragment.choose(ms)
    sm2=Chem.MolToSmiles(largest_mol)
    return sm2

def main1():
    df=file1
    try:
       df['Conc 1 Mean Op (Standardized)']=df['Conc 1 Mean Op (Standardized)'].fillna('=')
    except:
       messagebox.showinfo('Error','There is no column with heading Conc 1 Mean Op (Standardized), file can not be processed')
    inp=file2
    inp=inp.fillna('NA')
    l1=[i for i in inp.columns.tolist() if len(inp[i].unique())>1]
    d1=process(df,inp,l1)
    pp=d1.columns.tolist()
    df5=d1.dropna()
    cut=5
    mk='Conc 1 Mean Op (Standardized)' #with values like '=',">"
    ml='Conc 1 Mean (Standardized)' #with numerical activity values
    df6=df5[df5[ml]!='NA']
    try:
       eq=df6[df6[mk]=="="]
       #eq=eq[eq[ml]!='NA']
       eq['Active']=eq.apply(lambda x: 1 if x[ml]<=cut else 0, axis=1)
    except:
        print('No = in Conc 1 Mean Op (Standardized)')
        #continue
    try:
       ap=df6[df6[mk]=='~']
       ap['Active']=ap.apply(lambda x: 1 if x[ml]<=cut else 0, axis=1)
    except:
        print('No ~ in Conc 1 Mean Op (Standardized)')
        #continue
    try:
       gt=df6[df6[mk]=='>']
       gtc=gt[gt[ml]>=cut]
       gtc['Active']=0
    except:
        print('No > in Conc 1 Mean Op (Standardized)')
        #continue
    try:
       lt=df6[df6[mk]=='<']
       ltc=lt[lt[ml]<=cut]
       ltc['Active']=1
    except:
        print('No < in Conc 1 Mean Op (Standardized)')
        #continue
    ct2=pd.concat([eq,ap,gtc,ltc],axis=0)
    #ct2['Smiles']=ct2.apply(lambda x:castosmi(str(x['CAS Number'])), axis=1)
    #ct2.shape
    ct2.to_csv('FirstFile_forChecking.csv', index=False)

def main2():
    sm=file4
    print(sm.iloc[:,1:2].columns[0])
    ct2=file3
    inp2=file5
    ct3=pd.merge(ct2,sm,on='CAS Number', how='left')
    ct3=ct3.dropna()
    ct3['Can_smi']=ct3.apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x[sm.iloc[:,1:2].columns[0]])), axis=1)
    ct3['Can_smiPr']=ct3.apply(lambda x:processSmi(x['Can_smi']), axis=1)
    li=inp2.columns.tolist()
    li.append('Can_smiPr')
    duplicate = ct3[ct3.duplicated(li[1:], keep=False)]
    #duplicate = ct3.drop_duplicates( subset = li[1:], keep = 'last')
    #ct3.shape[0]-duplicate.shape[0]
    nodup=ct3.drop(duplicate.index.tolist(), axis=0)
    nodup['Nature']='No duplicate'
    k1,k2=[],[]
    ls=li[1:]
    for name,group in duplicate.groupby(li[1:]):
        if len(group['Active'].unique())>1:
           k1.append(group)
           dupconf=pd.concat(k1,axis=0)
        else:
           k2.append(group)
           dupwoc=pd.concat(k2,axis=0)
    ls=ls+['Active']
    #dupwoc[ls].drop_duplicates()
    dupconf2=dupconf[ls]
    #dupwoc2=dupwoc[ls]
    dupwoc3=dupwoc.drop_duplicates( subset = ls, keep = 'last').reset_index(drop = True)
    dupwoc3['Nature']='Duplicates'
    ndp=pd.concat([nodup,dupwoc3], axis=0)
    ndp.to_csv('NoDuplicates.csv', index=False)
    #ndp.dropna().shape
    dupconf.to_csv('DupConflict.csv', index=False)


form = tk.Tk()

form.title("Ecotox Data Curator")
form.geometry("650x200")

tab_parent = ttk.Notebook(form)

tab1 = ttk.Frame(tab_parent)
tab2 = ttk.Frame(tab_parent)



tab_parent.add(tab1, text="Initial data generation")
tab_parent.add(tab2, text="Final data generation")



# === WIDGETS FOR TAB ONE
#first
T1L1 = tk.Label(tab1, text="Select ECOTOX downloaded file",font=("Helvetica", 12))
T1L1.place(x=20,y=10)
T1E1 = tk.Entry(tab1,text='',width=40)
T1E1.place(x=260,y=13)
T1B1=tk.Button(tab1,text='Browse', command=data1,font=("Helvetica", 10))
T1B1.place(x=520,y=10)

T1L2 = tk.Label(tab1, text="Select first input file",font=("Helvetica", 12))
T1L2.place(x=113,y=50)
T1E2 = tk.Entry(tab1,text='',width=40)
T1E2.place(x=260,y=53)
T1B2=tk.Button(tab1,text='Browse', command=data2,font=("Helvetica", 10))
T1B2.place(x=520,y=50)

T1L3 = tk.Label(tab1, text="Activity cutoff in Conc 1 Mean (Standardized)",font=("Helvetica", 12),anchor=W, justify=LEFT)
T1L3.place(x=50,y=100)
T1E3 = tk.Entry(tab1)
T1E3.place(x=370,y=100)

B1=Button(tab1, text='Generate first data', command=main1,bg="orange",font=("Helvetica", 10),anchor=W, justify=LEFT)
B1.place(x=310,y=130)

# === WIDGETS FOR TAB TWO
T2L1 = tk.Label(tab2, text="Select modified ECOTOX file",font=("Helvetica", 12))
T2L1.place(x=20,y=10)
T2E1 = tk.Entry(tab2,text='',width=40)
T2E1.place(x=260,y=13)
T2B1=tk.Button(tab2,text='Browse', command=data3,font=("Helvetica", 10))
T2B1.place(x=520,y=10)

T2L2 = tk.Label(tab2, text="Select SMILES file",font=("Helvetica", 12))
T2L2.place(x=113,y=50)
T2E2 = tk.Entry(tab2,text='',width=40)
T2E2.place(x=260,y=53)
T2B2=tk.Button(tab2,text='Browse', command=data4,font=("Helvetica", 10))
T2B2.place(x=520,y=50)

T2L3 = tk.Label(tab2, text="Select second input file",font=("Helvetica", 12))
T2L3.place(x=93,y=90)
T2E3 = tk.Entry(tab2,text='',width=40)
T2E3.place(x=260,y=93)
T2B3=tk.Button(tab2,text='Browse', command=data5,font=("Helvetica", 10))
T2B3.place(x=520,y=90)


B2=Button(tab2, text='Generate final data', command=main2,bg="orange",font=("Helvetica", 10),anchor=W, justify=LEFT)
B2.place(x=310,y=130)



tab_parent.pack(expand=1, fill='both')

form.mainloop()
