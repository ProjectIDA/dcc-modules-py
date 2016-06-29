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

