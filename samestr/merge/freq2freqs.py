from os.path import isfile
import logging
import numpy as np

from samestr.utils.utilities import load_numpy_file

LOG = logging.getLogger(__name__)


def merge_freqs(freqs, freq, freq_name):
    """Appends two numpy arrays."""

    if freqs.shape[1] == freq.shape[1]:
        _freqs = np.append(freqs, freq, axis=0)
        if _freqs.shape[0] == freqs.shape[0] + freq.shape[0]:
            return _freqs, True
        else:
            LOG.error('Could not add array: %s' % freq_name)
    else:
        LOG.error('Dimensions of arrays do not match: %s (%s|%s)' %
                  (freq_name, freqs.shape[1], freq.shape[1]))

    return freqs, False


def freq2freqs(args):
    """Merges nucleotide frequencies for individual clades."""

    clade_freqs = None
    samples = []

    # skip if clade file exists
    output_file = '%s/%s' % (args['output_dir'], args['clade'])
    if isfile('%s.npz' % output_file):
        LOG.info('Skipping: %s. Output file existed.' % args['clade'])
        return True
    
    if len(args['input_files']) > 1:
        LOG.info('Merging %s inputs: %s' %
            (len(args['input_files']), args['clade']))

    # iterate samples containing clade
    for i, (sample, file_path) in enumerate(args['input_files']):

        multiple_samples = isinstance(sample, list)

        sample_clade_freqs = load_numpy_file(file_path)

        if clade_freqs is None:
            # final clade_freqs array will have dimensions (n_samples x clade_positions x 4)
            clade_freqs = np.zeros((len(args['input_files']), sample_clade_freqs.shape[1], sample_clade_freqs.shape[2]))
        else:
            if sample_clade_freqs.shape[1:] != clade_freqs.shape[1:]:
                LOG.error('Dimensions of arrays do not match: %s (%s|%s)' %
                  (sample, sample_clade_freqs.shape[1], clade_freqs.shape[1]))
                return False

        clade_freqs[i] = sample_clade_freqs

        if multiple_samples:
            samples += sample
        else:
            samples.append(sample)


        # # initialize clade freqs arrays
        # if clade_freqs is None:
        #     clade_freqs = load_numpy_file(file_path)
        #     if multiple_samples:
        #         samples += sample
        #     else:
        #         samples.append(sample)

        #     # if more than one sample, continue samples loop
        #     if len(args['input_files']) > 1:
        #         LOG.info('Merging %s inputs: %s' %
        #                  (len(args['input_files']), args['clade']))
        #         continue

        # if more than one sample, add sample freq to freqs
        # if len(args['input_files']) > 1:
        #     clade_freqs, success = merge_freqs(
        #         clade_freqs, load_numpy_file(file_path), sample)
        #     if success:
        #         if multiple_samples:
        #             samples += sample
        #         else:
        #             samples.append(sample)

    # per clade, save freqs and names to file
    np.savez_compressed(output_file + '.npz', clade_freqs, allow_pickle=True)

    with open(output_file + '.names.txt', 'w') as file:
        file.write('\n'.join(samples) + '\n')
