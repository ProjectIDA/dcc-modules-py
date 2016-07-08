from functools import reduce
from ida.signals.paz import PAZ

def pick(picklist, title=None, prompt=None, allow_quit_q=False, menu_on_error=False, err_message=None, indent_width=4):

    indent = indent_width*' '
    list_item_fmt = '{:>'+str(indent_width+1)+'}) {}'
    ndxlist = list(range(1, len(picklist)+1))

    quit = False
    display_list = True
    choice = 0
    choice_list = [str(inum) for inum in ndxlist]

    if not prompt:
        prompt = 'Enter number of selection: '
    prompt = indent + prompt

    while not quit and (choice not in choice_list):

        print('\n')
        if title:
            print(indent + title)
            print(indent + len(title)*'-', '\n')
        else:
            print()

        if display_list:
            for ndx, opt in enumerate(picklist):
                print(list_item_fmt.format(ndx + 1, opt))
        print()

        choice = input(prompt)

        if allow_quit_q and (choice in ['q', 'Q']):
            quit = True
            ndx = 0
            break
        elif choice not in choice_list:
            if err_message:
                print('ERROR\nERROR: ', err_message, '\nERROR')
            display_list = menu_on_error
        else:
            ndx = int(choice)
            break

    return quit, ndx-1

def pick2(picklistgroups, title=None, group_titles=None, prompt=None, multiple_choice=False,
          implicit_quit_q=False, menu_on_error=False, err_message=None, indent_width=4):

    # print('Pick list groups supplied:', picklistgroups)

    indent = indent_width*' '
    list_item_fmt = '{:>'+str(indent_width+1)+'}) {}'

    display_groups = []
    valid_ndx_tpls = {}
    user_choices = []
    user_choice_groups = {}  # dict to hold lists of index choices within each group; keyed by group index (int)

    choice_ndx = 1  # starting choice index for group lsits without user supplied choice keys
    for gndx, picklist in enumerate(picklistgroups):
        user_choice_groups[gndx] = []  # initialize user choices to empty ndx list for each group
        agroup = []
        for iindx, item in enumerate(picklist):
            if isinstance(item, str):
                # print('Index item:', item)
                choice_key = str(choice_ndx)
                choice_ndx += 1
                agroup.append((choice_key, item))
                valid_ndx_tpls[choice_key] = (gndx, iindx)
            elif isinstance(item, tuple):
                agroup.append((item[0], item[1]))
                valid_ndx_tpls[item[0].upper()] = (gndx, iindx)
            else:
                print('uh oh:', iindx, item)

        display_groups.append(agroup)

    if not prompt:
        prompt = 'Enter selection: '
    prompt = indent + prompt

    user_quit = False
    is_valid = False
    display_list = True

    while (not user_quit) and (not is_valid):

        print('\n')
        if title:
            print(indent + title)
            print(indent + len(title)*'-')
        else:
            print()

        if display_list:

            for gndx, glist in enumerate(display_groups):
                if group_titles:
                    print()
                    if len(group_titles[gndx].strip()) > 0:
                        print(indent + group_titles[gndx])
                if len(glist) > 0:
                    for choice in glist:
                        print(list_item_fmt.format(choice[0], choice[1]))
                else:
                    print(indent + '(none)')

        print()

        choice = input(prompt).upper().strip()
        user_choices = [chc.strip() for chc in choice.split(',')]
        user_quit = ('Q' in user_choices) and implicit_quit_q

        # print('pick2: choice:', choice)
        # print('pick2: user_choices:', user_choices)
        if multiple_choice:
            is_valid = reduce(lambda x, y: x and (y in valid_ndx_tpls.keys()), user_choices, True)
        else:
            is_valid = (len(user_choices) > 0) and (user_choices[0] in valid_ndx_tpls.keys())

        if is_valid:
            user_quit = False   # 'Q' must be a valid option, so override possible implicit_quit
            for chc in user_choices:
                user_choice_groups[valid_ndx_tpls[chc][0]].append(valid_ndx_tpls[chc][1])
            break
        elif (not user_quit):
            if err_message:
                print('\nERROR\nERROR: ', err_message, '\nERROR')
            display_list = menu_on_error
        else:  # user_quit must == True
            break

    return (not user_quit), user_choices, user_choice_groups


def select_perturb_map(paz):

    if not isinstance(paz, PAZ):
        raise TypeError('paz must be a populated PAZ object')

    paz_fit_lf = paz.make_partial2(norm_freq=1.0, partial_mode=paz.PARTIAL_FITTING_LF)
    paz_fit_hf = paz.make_partial2(norm_freq=1.0, partial_mode=paz.PARTIAL_FITTING_HF)

    poles_pert_def, zeros_pert_def = paz.perturb_defaults()
    defchoice = [('D', 'Use Defaults (indicated by "<=="')]

    # make list of LF p/z to user perturbing choices
    # Do LOW FREQ FIRST
    zero_pert_choices = []
    pole_pert_choices = []
    for ndx, val in enumerate(paz_fit_lf.zeros()):
        if ndx in zeros_pert_def[0]:
            zero_pert_choices.append(str(val) + ' <==')
        else:
            zero_pert_choices.append(str(val))
    for ndx, val in enumerate(paz_fit_lf.poles()):
        if ndx in poles_pert_def[0]:
            pole_pert_choices.append(str(val) + ' <==')
        else:
            pole_pert_choices.append(str(val))

    pert_choices = [defchoice, zero_pert_choices, pole_pert_choices]
    success, choices, pert_choice_groups = pick2(pert_choices, 'Select LOW Freq zeros & poles to perturb',
                                                 prompt='Enter selection (or "q" to quit): ',
                                                 group_titles=['',
                                                               'LOW Freq Zeros',
                                                               'LOW Freq Poles'],
                                                 multiple_choice=True,
                                                 implicit_quit_q=True, menu_on_error=True)

    # print(success, pert_choice_groups)
    if not success:
        return False, None, None

    if choices[0].upper() == 'D':  # using defaults
        lf_map = (poles_pert_def[0], zeros_pert_def[0])  # beware, put poles then zeros in this map tuple
    else:
        lf_map = (pert_choice_groups[2], pert_choice_groups[1])


    # NOW HIGH FREQ
    zero_pert_choices = []
    pole_pert_choices = []
    for ndx, val in enumerate(paz_fit_hf.zeros()):
        if ndx in zeros_pert_def[1]:
            zero_pert_choices.append(str(val) + ' <==')
        else:
            zero_pert_choices.append(str(val))
    for ndx, val in enumerate(paz_fit_hf.poles()):
        if ndx in poles_pert_def[1]:
            pole_pert_choices.append(str(val) + ' <==')
        else:
            pole_pert_choices.append(str(val))

    pert_choices = [defchoice, zero_pert_choices, pole_pert_choices]
    success, choices, pert_choice_groups = pick2(pert_choices, 'Select HIGH Freq zeros & poles to perturb',
                                                 prompt='Enter selection (or "q" to quit): ',
                                                 group_titles=['',
                                                               'HIGH Freq Zeros',
                                                               'HIGH Freq Poles'],
                                                 multiple_choice=True,
                                                 implicit_quit_q=True, menu_on_error=True)

    # print(success, pert_choice_groups)

    if not success:
        return False, None, None

    if choices[0].upper() == 'D':  # using defaults
        hf_map = (poles_pert_def[1], zeros_pert_def[1])  # beware, put poles then zeros in this map tuple
    else:
        hf_map = (pert_choice_groups[2], pert_choice_groups[1])


    return success, lf_map, hf_map  # each in (poles, zeros) order
