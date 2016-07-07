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

        print()
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
                print('\nERROR\nERROR: ', err_message, '\nERROR')
            display_list = menu_on_error
        else:
            ndx = int(choice)
            break

    return quit, ndx-1

def pick2(picklistgroups, title=None, group_titles=None, prompt=None, allow_quit_q=False, menu_on_error=False, err_message=None, indent_width=4):

    # print('Pick list groups supplied:', picklistgroups)

    indent = indent_width*' '
    list_item_fmt = '{:>'+str(indent_width+1)+'}) {}'

    choice_cnt = reduce(lambda x, y: x + len(y), picklistgroups, 0)
    ndxlist = list(range(1, choice_cnt+1))

    groups = []
    valid_choices = {}
    choice_ndx = 1
    for gndx, picklist in enumerate(picklistgroups):
        agroup = []
        for iindx, item in enumerate(picklist):
            if isinstance(item, str):
                # print('Index item:', item)
                choice_key = str(choice_ndx)
                choice_ndx += 1
                agroup.append((choice_key, item))
                valid_choices[choice_key] = (gndx, iindx)
            elif isinstance(item, tuple):
                print('Tuple item:', picklist)
                agroup.append((item[0], item[1]))
                valid_choices[item[0].upper()] = (gndx, iindx)
            else:
                print('uh oh:', iindx, item)

        print('this group:', agroup)
        groups.append(agroup)
        print('all groups:', groups)

    quit = False
    display_list = True
    choice = None

    print('Valid choices:', valid_choices)

    if not prompt:
        prompt = 'Enter number of selection: '
    prompt = indent + prompt
            
    while not quit and (choice not in valid_choices):

        print()
        if title:
            print(indent + title)
            print(indent + len(title)*'-')
        else:
            print()

        if display_list:

            for gndx, glist in enumerate(groups):
                if group_titles:
                    print()
                    if len(group_titles[gndx].strip()) > 0:
                        print(indent + group_titles[gndx])
                for choice in glist:
                    print(list_item_fmt.format(choice[0], choice[1]))

        print()

        choice = input(prompt).upper()

        if choice in valid_choices:
            quit = False
            restpl = valid_choices[choice]
            itemtxt = picklistgroups[restpl[0]][restpl[1]]
            break
        if allow_quit_q and (choice in ['q', 'Q']):
            choice = 'q'
            quit = True
            restpl = None
            itemtxt = None
            break
        else:
            if err_message:
                print('\nERROR\nERROR: ', err_message, '\nERROR')
            display_list = menu_on_error

    return not quit, choice, itemtxt, restpl

