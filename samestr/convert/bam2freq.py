from os.path import isfile
from os import remove
import logging

from glob import glob
from samestr.utils import ooSubprocess

from samestr.convert.kpileup import pileup
from samestr.convert.kp2np import kp2np

LOG = logging.getLogger(__name__)


def bam2freq(arg):
    with open(arg["kp"], "wt") as kpileup_out:
        #  sample_id, bam_file, gene_file, min_bq, min_mq, min_depth, outstream=kpileup_out
        position_stream = pileup(
            arg['bname'],
            arg['sorted_bam'],
            arg['gene_file'],
            arg['min_base_qual'],
            arg['min_aln_qual'],
            arg['min_vcov'],
            outstream=kpileup_out,
        )
        kp2np(
            position_stream,
            arg['contig_map'],
            arg['bname'],
            arg['gene_file'],
            arg['np'],
        )
        if not arg['keep_tmp_files']:
            remove(arg['gene_file'])
            remove(arg['contig_map'])
            remove(arg['kp'])
            remove(arg['sorted_bam'])
            remove(arg['sorted_bam'] + '.bai')

        LOG.debug('Finished: %s' % arg['np'])
        
    return arg

def bam2freq_old(arg):
    """
        converts mapping files
        from binary sequence alignment (bam) format
        to SNV profiles (numpy arrays).
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])
    if not arg:
        LOG.error('Empty input')
        return None

    if not isfile(arg['kp']):
        LOG.debug('Piling: %s' % arg['sorted_bam'])
        oosp.ex(
            'kpileup.py',
            args=[
                arg['bname'],
                arg['sorted_bam'],
                arg['gene_file'],
                str(arg['min_base_qual']),
                str(arg['min_aln_qual']),
                str(arg['min_vcov'])
            ],
            out_fn=arg['kp'],
            verbose=False)

    if not glob('%s/*.np.cPickle' % arg['np']) or \
       not glob('%s/*.npy' % arg['np']) or \
       not glob('%s/*.npz' % arg['np']) or \
       not glob('%s/*.npy.gz' % arg['np']):
        LOG.debug('Converting: %s' % arg['kp'])
        oosp.ex('kp2np.py',
                args=[
                    '--kp', arg['kp'], '--map', arg['contig_map'], '--sample',
                    arg['bname'], '--gene-file', arg['gene_file'],
                    '--output-dir', arg['np']
                ],
                verbose=False)
    else:
        LOG.warning('Skipping: %s. '
                    'Directory already contained np files: %s' %
                    (arg['bname'], arg['np']))
        
    # clean up
    if not arg['keep_tmp_files']:
        remove(arg['gene_file'])
        remove(arg['contig_map'])
        remove(arg['kp'])
        remove(arg['sorted_bam'])
        remove(arg['sorted_bam'] + '.bai')

    LOG.debug('Finished: %s' % arg['np'])
    return arg
