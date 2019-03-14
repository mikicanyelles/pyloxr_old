from tkinter import *
#from old_modules import *
import old_modules as om

root=Tk()
root.title('pyLOXr')
#root.iconbitmap('favicon.ico')
#root.config(bg="black")

#main_frame=Frame(root)
#main_frame.pack(fill="both", expand=True)
#main_frame.config(width=150, height=180)

#md_sum=lambda:om.summarize(top, dyn)
def md_sum():
    print('holiii')

top='3rde.prmtop'; dyn='5_prod_09.nc'
but1=Button(root, text=" MD Summariser ", command=md_sum, fg='red')#(top,dyn))#, command=old_modules.summarize())
but1.grid(row=0, column=0, padx=20, pady=10)
#but1.set(old_modules.summarize_test())

butExit=Button(root, text=" Exit ", command=quit)
butExit.grid(row=3, column=1, padx=10, pady=10)
#entry1=Entry(main_frame)
#entry1.grid(row=1, column=0)


root.mainloop()
