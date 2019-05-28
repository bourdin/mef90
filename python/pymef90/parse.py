import yaml
import argparse
def parse(parser,key=None):
    args = parser.parse_args()

    if key:
        arg_dict = args.__dict__
        if arg_dict[key]:
            print('WARNING: reading options from YAML file. Duplicate and positional command line options will be ignored. \n')
            data = yaml.load(arg_dict[key],Loader=yaml.FullLoader)
            # replace the file object with its name
            arg_dict[key] = arg_dict[key].name
            # optional and positional arguments
            for group in parser._action_groups:
                if group.title in ['optional arguments', 'positional arguments']:
                    for key, value in data.items():
                        if not isinstance(value, dict):
                            arg_dict[key] = value

            # groups
            for group in parser._action_groups:
                if group.title in data.keys():
                    for key, value in data[group.title].items():
                        if isinstance(value, list):
                            for v in value:
                                arg_dict[key].append(v)
                        else:
                            arg_dict[key] = value
            # replace the file object with its name
            #arg_dict[key] = arg_dict[key].name
    return args

def parseGroup(parser,args,groupName): 
    for group in parser._action_groups:
        if group.title == groupName:
            return argparse.Namespace(**{a.dest:getattr(args,a.dest,None) for a in group._group_actions})
    return argparse.Namespace()
