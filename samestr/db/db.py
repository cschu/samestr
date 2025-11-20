import logging
import pathlib

from os.path import isfile

import pandas as pd

from samestr.db import generate_db
from samestr.db.db_checks import check_args, check_db_integrity
from samestr.utils.file_mapping import get_uniform_extension
from samestr.utils.utilities import read_json, read_tsv


LOG = logging.getLogger(__name__)

def load_samestr_db_manifest(db_path, samestr_cmd):
    """
    read the samestr database files:
    - db_manifest.json
    - db_clades.json
    - db_taxonomy.tsv
    """
    db_path = pathlib.Path(db_path)

    db_manifest_fn = db_path / 'db_manifest.json'
    db_clades_fn = db_path / 'db_clades.json'
    db_taxonomy_fn = db_path / 'db_taxonomy.tsv'

    db_existed = db_manifest_fn.is_file()

    if db_existed:
        db_manifest = read_json(db_manifest_fn)
        db_clades = read_json(db_clades_fn)
        db_taxonomy = {
            'fpath': db_taxonomy_fn,
            'records': read_tsv(db_taxonomy_fn)
        }
    else:
        db_manifest = {
            'fpath': db_manifest_fn,
            'database': {},
            'records': {
                'total_n_files': 0,
                'total_n_clades': 0,
                'total_n_markers': 0,
                'total_n_positions': 0,
            }
        }
        db_clades = {
            'fpath': db_clades_fn,
            'records': {}
        }
        db_taxonomy = {
            'fpath': db_taxonomy_fn,
            'records': pd.DataFrame(
                columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'clade']
            )
        }
        
    return {
        "db_existed": db_existed,
        "db_manifest": db_manifest,
        "db_clades": db_clades,
        "db_taxonomy": db_taxonomy,
        "samestr_cmd": samestr_cmd,
    }


def load_db(input_args, samestr_cmd):
    ss_db = load_samestr_db_manifest(
        input_args['marker_dir'],
        samestr_cmd
    )

    if ss_db['db_existed']:
        db_manifest = ss_db['db_manifest']
        LOG.info(
            'Processing data with existing SameStr database based on %s, %s' % \
                (db_manifest['database']['name'], db_manifest['database']['version'])
        )
        LOG.info(
            'Database contains %s files, %s clades, %s markers, totalling %s positions.' % \
                (
                    db_manifest['records']['total_n_files'],
                    db_manifest['records']['total_n_clades'],
                    db_manifest['records']['total_n_markers'],
                    db_manifest['records']['total_n_positions']
                )
        )
    elif not ss_db['db_existed'] and input_args['command'] != 'db':
        LOG.error(
            'Aborting. Could not parse database at %s' % \
                input_args['output_dir']
            )
        exit(1)

    return ss_db


def db_main(input_args, samestr_cmd, db_extensions):
    check_args(input_args)

    load_db(input_args, samestr_cmd)

    #### check db integrity if requested and db existed
    if input_args.get("db_check"):
        check_db_integrity(input_args)

    input_args['input_extension'] = get_uniform_extension(
        [input_args['markers_fasta']],
        db_extensions,
    )

    output_dir = pathlib.Path(input_args['output_dir']) / "db_markers"
    output_dir.mkdir(exist_ok=True, parents=True,)

    # expand and generate db from markers/info files
    generate_db(input_args)
