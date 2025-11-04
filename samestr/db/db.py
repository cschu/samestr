import pathlib

from samestr.db import generate_db
from samestr.db.db_checks import check_args, check_db_integrity
from samestr.db.load_db import load_db
from samestr.utils.file_mapping import get_uniform_extension


accepted_extensions = ['.fasta', '.fa', '.fna', '.fasta.bz2', '.fa.bz2', '.fna.bz2']

def db_main(input_args, samestr_cmd):
	check_args(input_args)
	
	load_db(input_args, samestr_cmd)

	#### check db integrity if requested and db existed
	if input_args.get("db_check"):
		check_db_integrity(input_args)

	input_args['input_extension'] = get_uniform_extension(
		[input_args['markers_fasta']],
		accepted_extensions,
	)

	output_dir = pathlib.Path(input_args['output_dir']) / "db_markers"
	output_dir.mkdir(exist_ok=True, parents=True,)

    # expand and generate db from markers/info files
	generate_db(input_args)