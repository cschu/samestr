import logging
import os
import pathlib

from os.path import isfile, isdir

from samestr.utils.file_mapping import clade_path


LOG = logging.getLogger(__name__)


def check_args(input_args):
    if not input_args['db_check']:
        for o in ['db_version', 'markers_fasta', 'markers_info']:
            if o not in input_args:
                LOG.error('The following argument is required: %s' % o)
                exit(1)
        input_args['marker_dir'] = input_args['output_dir']


def check_db_integrity(input_args):
    ss_db = input_args['samestr_db']
    if not ss_db['db_existed']:
        LOG.error('Cannot conduct a SameStr database check since a database was not provided.')
        exit(1)
    else:
        t_missing = []
        t_found = 0
        LOG.info('Checking integrity of the SameStr database.')
        expected_file_bases = ['%s/%s' % (input_args['marker_dir'], clade_path(c, filebase = True)) for c in ss_db['db_clades']['records'].keys()]
        for e in  expected_file_bases:
            suffix = ['.contig_map.txt.gz', '.gene_file.txt.gz', '.markers.fa.gz', '.positions.txt.gz']
            for s in suffix:
                t = e + s
                if not isfile(t):
                    t_missing.append(t)
                else:
                    t_found += 1

        if t_found == ss_db['db_manifest']['records']['total_n_files']:
            LOG.info('Number of files in the SameStr database matches the manifest.')
        if len(t_missing) > 0:
            outfn = os.path.join(input_args['output_dir'], 'db_check.tsv')
            LOG.warning('Files missing from the SameStr database: %s. Saving to %s' % (len(t_missing), outfn))
            pathlib.Path(input_args['output_dir']).mkdir(exist_ok=True, parents=True,)
            with open(outfn, 'w') as f:
                print(*t_missing, file=f, sep="\n",)
                
        exit(0)