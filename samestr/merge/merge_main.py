from os.path import basename, isfile

from samestr.merge.freq2freqs import freq2freqs
from samestr.utils.ooSubprocess import parallelize_async


def merge_main(input_args, input_clade):
	clade_file_dict = {}

	# group input files by clade
	for file_path in input_args['input_files']:
		fn = basename(file_path)
		clade = fn.split('.')[0]

		# append list if file has been merged before
		sample_names = file_path.replace(
			input_args['input_extension'], '.names.txt'
		)
		if isfile(sample_names):
			with open(sample_names, 'r', encoding="utf8") as s_n:
				sample = s_n.read().strip().split('\n')
		else:
			sample = fn.split(input_args['input_extension'])[0].replace(f'{clade}.', '')

		# skip if not in selected clade
		# if input_clade is not None and clade not in input_clade:
		#     continue
		if input_clade is None or clade in input_clade:
			clade_file_dict.setdefault(clade, []).append([sample, file_path])

	cmd_args = [
		{
			"clade": clade,
			"input_files": files,
			"output_dir": input_args["output_dir"],
		}
		for clade, files in clade_file_dict.items()
	]        

	cmd_args = parallelize_async(freq2freqs, cmd_args, input_args['nprocs'])
