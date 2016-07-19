from functools import reduce


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
