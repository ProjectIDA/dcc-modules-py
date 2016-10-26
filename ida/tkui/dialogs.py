from tkinter import font, Toplevel, Frame, Label, Button, NW, LEFT

class AboutDialog(Toplevel):
    # A dialog box to get reistration number
    # This is not needed unless the car was speeding
    def __init__(self, master, about_head='About', about_text='About', *args, **kwargs):
        Toplevel.__init__(self, master, *args, **kwargs)

        self.resizable(width=False, height=False)

        # position window over parent
        # self.geometry("350x800+%d+%d" % (master.winfo_rootx(),
        #               master.winfo_rooty()))

        self.title(about_head.replace('\n', ' '))
        main_frm = Frame(self, padx=10, pady=10)
        main_frm.grid()
        main_frm.rowconfigure(0, weight=1)
        main_frm.rowconfigure(1, weight=1)
        main_frm.rowconfigure(2, weight=1)
        hdlbl = Label(main_frm, text=about_head, anchor=NW, justify=LEFT)
        hdlbl.grid(row=0, column=0, sticky='w')
        hdfont = font.Font(font=hdlbl['font'])
        hdfont.config(weight='bold')
        hdlbl.config(font=hdfont)
        Label(main_frm, text=about_text, anchor=NW, justify=LEFT).grid(row=1, column=0, sticky='w')
        Button(main_frm, text="OK", command=self.on_ok).grid(row=2, column=0)

        # make dialog modal
        self.transient(master)
        self.focus_set()
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self.on_ok)

    def on_ok(self):
        # called when ok button is pressed
        self.destroy()
