from samestr.convert import concatenate_gene_files, bam2freq, sam2bam
from samestr.utils.file_mapping import set_output_structure, spread_args_by_input_files
from samestr.utils.ooSubprocess import parallelize_async
from samestr.utils.utilities import list_group_by_basename


def convert_main(input_args):
	# TODO: offer/fix dir-based input
        # input_args['input_files'], input_args['input_extension'] = collect_input_files(
        #     input_args['input_dir'],
        #     accepted_extensions,
        #     recursive=input_args["recursive_input"]
        # )
        
        input_args['input_files'] = list_group_by_basename(
            input_args['input_files'],
            cut_name_endings=[input_args['input_extension']])

        # Spread args over files/file-pairs, Set output names and check/make their dir
        input_args['input_sequence_type'] = 'single'
        cmd_args = spread_args_by_input_files(input_args)
        cmd_args = set_output_structure(cmd_args)

        # Run: convert
        cmd_args = parallelize_async(concatenate_gene_files, cmd_args,
                                 input_args['nprocs'])
        
        cmd_args = parallelize_async(sam2bam, cmd_args, input_args['nprocs'])        
        cmd_args = parallelize_async(bam2freq, cmd_args, input_args['nprocs'])
