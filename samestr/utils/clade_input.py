import logging

from os.path import basename

LOG = logging.getLogger(__name__)

def manage_clade_input(input_args, input_clade):
    # For each species: group freqs with resp names
    freqs = {
        basename(fn).split(input_args['input_extension'])[0]: [fn]
        for fn in input_args['input_files']
    }

    for name in input_args['input_names']:
        clade = basename(name).split('.names.txt')[0]
        if clade not in freqs:
            LOG.warning('Skipping %s. Found name file '
                        'but no SNV profile.', clade)
        else:
            freqs[clade].append(name)

    for clade, freq in list(freqs.items()):
        if len(freq) != 2:
            LOG.warning('Skipping [%s]. Found SNV profile '
                        'but no name file: %s' % (clade, ', '.join(freq)))
            freqs.pop(clade)
        elif input_clade is not None and clade not in input_clade:
            freqs.pop(clade)

    # attach sample selection file
    if 'samples_select' in input_args:
        for i_s in input_args['samples_select']:
            clade = basename(i_s).split('.select.txt')[0]
            if clade not in freqs:
                LOG.warning('Skipping %s. Found selection file '
                            'but no SNV profile.', clade)
            else:
                freqs[clade].append(i_s)

    for clade, freq in list(freqs.items()):
        if len(freq) > 3:
            LOG.warning('Skipping %s. Found more than one '
                        'sample selection file: %s' % (clade, ', '.join(freq)))
            freqs.pop(clade)
        elif len(freq) == 2:
            freq += [None]
        elif len(freq) != 3:
            LOG.error('Unexpected number of files (%s): %s.', len(freq), clade)
            exit(0)  

    # Spread args over clades
    cmd_args = []
    ignore_args = {'input_files', 'input_names', 'input_select', 'clade'}

    arg_template_d = {
        arg: value
        for arg, value in input_args.items()
        if arg not in ignore_args
    }

    for clade, (input_file, input_name, input_select) in freqs.items():
        
        args_d = {
            "input_file": input_file,
            "input_name": input_name,
            "input_select": input_select,
            "clade": clade,
        }

        args_d.update(arg_template_d)

        cmd_args.append(args_d)
        
    return cmd_args, freqs