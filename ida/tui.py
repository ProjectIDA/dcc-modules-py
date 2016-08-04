#######################################################################################################################
# Copyright (C) 2016  Regents of the University of California
#
# This is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License (GNU GPL) as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE.TXT in the root of the source code repository.
# Additionally, it can be found at http://www.gnu.org/licenses/.
#
# NOTES: Per GNU GPLv3 terms:
#   * This notice must be kept in this source file
#   * Changes to the source must be clearly noted with date & time of change
#
# If you use this software in a product, an explicit acknowledgment in the product documentation of the contribution
# by Project IDA, Institute of Geophysics and Planetary Physics, UCSD would be appreciated but is not required.
#######################################################################################################################
from functools import reduce
from enum import Enum

from fabulous.color import bold, underline, blue, red

class PickResult(Enum):
    collect_noop = 0
    collect_ok = 1
    collect_back = 2
    collect_error = 3
    collect_quit = 4


def pick(picklist, title=None, prompt=None, allow_quit_q=False, menu_on_error=False, err_message=None, indent_width=4):

    indent = (indent_width+1)*' '
    list_item_fmt = indent + '{} {}'
    ndxlist = list(range(1, len(picklist)+1))

    quit = False
    display_list = True
    choice = -1
    choice_list = [str(inum) for inum in ndxlist]

    if not prompt:
        prompt = 'Enter number of selection: '
    prompt = indent + prompt

    while not quit and (choice not in choice_list):

        print('\n')
        if title:
            print(indent + underline(bold(blue(title))))
            # print(indent + blue(len(title)*'-'), '\n')

        print()
        if display_list:
            for ndx, opt in enumerate(picklist):
                print(list_item_fmt.format(str(blue(bold(str(ndx + 1))+')')), str(bold(opt))))
        print()

        choice = input(bold(blue(prompt)))

        if allow_quit_q and (choice in ['q', 'Q']):
            quit = True
            ndx = 0
            break
        elif choice not in choice_list:
            if err_message:
                print('\n' + indent + red(bold(err_message)))
            display_list = menu_on_error
        else:
            ndx = int(choice)
            break

    return not quit, ndx-1


def pick2(picklistgroups, title=None, group_titles=None, prompt=None, multiple_choice=False,
          implicit_quit_q=False, implicit_back_b=False, menu_on_error=False, err_message=None, indent_width=4):

    indent = (indent_width+1)*' '
    list_item_fmt = indent + '{} {}'
    no_entry_message = 'No selection was entered. Please try again.'

    display_groups = []
    valid_ndx_tpls = {}
    user_choices = []
    user_choice_groups = []  # list to hold lists of index choices within each group; same seq as suuplied picklistgroup
    user_choice_tuples = []  # list to hold tuples of (grpndx, iindx) for each choice made by user;

    choice_ndx = 1  # starting choice index for group lsits without user supplied choice keys
    for gndx, picklist in enumerate(picklistgroups):
        # user_choice_groups[gndx] = []  # initialize user choices to empty ndx list for each group
        user_choice_groups.append([])  # initialize user choices to empty ndx list for each group
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
                print('Error processing picklist:', iindx, item)

        display_groups.append(agroup)

    if not prompt:
        prompt = 'Enter selection: '
    prompt = indent + prompt

    invalid = True
    display_list = True

    result = PickResult.collect_noop
    while (result not in [PickResult.collect_ok, PickResult.collect_back, PickResult.collect_quit]) and (invalid):

        print('\n')
        if title:
            print(indent + underline(bold(blue(title))))

        print()
        if display_list:
            for gndx, glist in enumerate(display_groups):
                if group_titles:
                    if gndx > 0:    # add blank line after previous group of choices
                        print()
                    if len(group_titles[gndx].strip()) > 0:
                        print(indent + underline(group_titles[gndx]))
                if len(glist) > 0:
                    for choice in glist:
                        # print(list_item_fmt.format(choice[0], choice[1]))
                        print(list_item_fmt.format(str(blue(bold(choice[0])) + ')'), str(bold(choice[1]))))
                else:
                    print(indent + '(none)')
        print()

        choice = input(bold(blue(prompt))).upper().strip()
        user_choices = [chc.strip() for chc in choice.split(',')]

        if user_choices and ('Q' == user_choices[0]) and implicit_quit_q:
            result = PickResult.collect_quit
        elif user_choices and ('B' in user_choices) and implicit_back_b:
            result = PickResult.collect_back
        elif not user_choices:
            print('\n' + indent + red(bold(no_entry_message)))
        else:
            if multiple_choice:
                invalid = not reduce(lambda x, y: x and (y in valid_ndx_tpls.keys()), user_choices, True)
            else:
                invalid = ((len(user_choices) != 1) or (user_choices[0] not in valid_ndx_tpls.keys()))

            if not invalid:
                result = PickResult.collect_ok
                for chc in user_choices:
                    user_choice_groups[valid_ndx_tpls[chc][0]].append(valid_ndx_tpls[chc][1])
                    user_choice_tuples.append(valid_ndx_tpls[chc])

            else:
                if err_message:
                    print('\n' + indent + red(bold(err_message)))
                display_list = menu_on_error

    return result, user_choices, user_choice_groups, user_choice_tuples

def input_yn(prompt, err_message=None, indent_width=4, default=None):

    indent = (indent_width+1)*' '
    if default:
        prompt = prompt + ' [{}]'.format(default)

    answered = False
    while not answered:

        answer = input(bold(blue(prompt))).upper().strip()
        if default and answer == '':
            answer = default
        answered = answer.upper() in ['Y', 'N']

        if not answered and err_message:
            print('\n' + indent + red(bold(err_message)) + '\n')

    return answer
