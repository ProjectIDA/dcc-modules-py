import re
from datetime import datetime
from tkinter import Entry

BGERR = '#FCC'
BGOK = '#CFC'
BGGRAY = '#CCC'

class SeedCodeEntry(Entry):
    """Used for entry single word text items like Network, Station"""

    def __init__(self, master=None, value='', charset_cre=None, **kw):
        Entry.__init__(self, master, **kw)
        self.insert('end', value)
        self._charset_cre = charset_cre

        vcmd = (self.register(self._validate), '%s', '%P')
        self.configure(validate='all', vcmd=vcmd)
        self.bind('<FocusOut>', self._focus_out)
        self._validate(value, value)


    def _focus_out(self, event):
        value = Entry.get(self)
        if value:
            newvalue = value.upper()
            self.delete(0, 'end')
            self.insert('end', newvalue)

    def _validate(self, old, new):
        if self._charset_cre:
            charset_ok = self._charset_cre.fullmatch(new) != None
        else:
            charset_ok = True

        return charset_ok or (new == '')


class ChanListEntry(Entry):

    chans_fil_re = '[a-zA-Z12, ]*'
    # chans_val_re = '(.hz( *, *| +).h1( *, *| +).h2)|(.hz( *, *| +).hn( *, *| +).he)'
    chans_fix_re = [r'( *, *| +)', ',', re.IGNORECASE]

    def __init__(self, master=None, value='', **kw):
        Entry.__init__(self, master, **kw)
        self.insert('end', value)
        self._valid = True
        vcmd = (self.register(self._validate), '%s', '%P')
        self.configure(validate='all', vcmd=vcmd)
        self.bind('<FocusOut>', self._focus_out)

        self._validate(value, value)

    def _focus_out(self, event):
        value = event.widget.get()
        if value and self._valid:
            value = value.strip().upper()
            newval = re.sub(self.chans_fix_re[0], self.chans_fix_re[1], value, flags=self.chans_fix_re[2])
            self.delete(0, 'end')
            self.insert('end', newval)

    def _validate(self, old, new):
        # self._valid = not (re.fullmatch(self.chans_val_re, new, flags=re.IGNORECASE) == None)
        filter = (re.fullmatch(self.chans_fil_re, new, flags=re.IGNORECASE) == None)
        # don't change if filtering'
        # if not filter:
        #     if self._valid:
        #         self.config(bg='#CFC')
        #     else:
        #         self.config(bg='#FCC')
        #
        return not filter


class IntEntry(Entry):

    _messages = {
        'en': {
            'below_min': 'Value [{}] is less than minimum permitted value [{}].',
            'above_max': 'Value [{}] is greater than maximum permitted value [{}].',
        }
    }

    def __init__(self, master=None, value=0, **kw):
        Entry.__init__(self, master, **kw)
        value = str(value)
        self.insert('end', value)
        # self._valid = True
        # self._minval = minvalue
        # self._maxval = maxvalue
        # self._errmsgs = []
        vcmd = (self.register(self._validate), '%s', '%P')
        self.configure(validate='all', vcmd=vcmd)

        self._validate(value, value)

    def _validate(self, old, new):
        new_valid = (new == '') or new.isdigit()
        # inrange = new_valid and (new != '')
        # self._errmsgs = []
        #
        # # if valid int
        # if new_valid:
        #     self._valid = new_valid
        #     curval = new
        # else:
        #     # if not valid int use previous value, re-assess validity adn range with old value
        #     self._valid = old.isdigit()
        #     curval = old
        #
        # # if valid int, whether prev or new value, check range.
        # if self._valid and (curval != ''):
        #     if self._minval:
        #         if not int(curval) >= self._minval:
        #             self._errmsgs.append(self._messages['en']['below_min'].format(curval, self._minval))
        #     if self._maxval:
        #         if not int(curval) <= self._maxval:
        #             self._errmsgs.append(self._messages['en']['above_max'].format(curval, self._maxval))
        #
        # if (not inrange) or (not self._valid) or (curval == ''):
        #     self.config(bg='#FCC')
        # else:
        #     self.config(bg='#CFC')
        #
        return new_valid # or (new == '')


class FloatEntry(IntEntry):

    float_val_re = '[0-9]*(\.){0,1}[0-9]*'

    def _validate(self, old, new):

        if new == '':
            new_valid = True
        else:
            try:
                float(new)
            except:
                new_valid = False
            else:
                new_valid = True

        # if new == '.': new = '0.'
        # new_valid = (re.fullmatch(self.float_val_re, new) != None)
        # inrange = True
        #
        # self._errmsgs = []
        # # if valid Float
        # if new_valid:
        #     self._valid = new_valid
        #     curval = new
        # else:
        #     # if not valid float use previous value, re-assess validity adn range with old value
        #     if old == '.': old = '0.'
        #     self._valid = not (re.fullmatch(self.float_val_re, old) == None)
        #     curval = old
        #
        # # if valid float, whether prev or new value, check range.
        # if self._valid and (curval != ''):
        #     if self._minval:
        #         if not float(curval) >= self._minval:
        #             self._errmsgs.append(self._messages['en']['below_min'].format(curval, self._minval))
        #     if self._maxval:
        #         if not float(curval) <= self._maxval:
        #             self._errmsgs.append(self._messages['en']['above_max'].format(curval, self._maxval))
        # flag error in UI of not valid float, not in range, or blank
        # (regardless of new or old value used)
        # if (not inrange) or (not self._valid) or (curval == ''):
        #     self.config(bg='#FCC')
        # else:
        #     self.config(bg='#CFC')
        #
        return new_valid


class IsoTimeEntry(Entry):

    isotime_fil_re = '[0-9:T\- ]*'

    def __init__(self, master=None, value='', **kw):
        Entry.__init__(self, master, **kw)
        self.insert('end', value)
        self._valid = True
        vcmd = (self.register(self._validate), '%s', '%P')
        self.configure(validate='all', vcmd=vcmd)
        # self.bind('<FocusOut>', self._focus_out)

        self._validate(value, value)

    def valid_isotime_string(self, timestr):
        try:
            dt = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')
        except:
            try:
                dt = datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S')
            except:
                valid = False
            else:
                valid = True
        else:
            valid = True

        return valid


    def _validate(self, old, new):
        self._valid = self.valid_isotime_string(new) and (len(new) == 19)
        filter = (re.fullmatch(self.isotime_fil_re, new) == None)

        self._errmsgs = []

        if (not self._valid):
            self.config(bg='#FCC')
        else:
            self.config(bg='#CFC')

        return (not filter) or (new == '')

