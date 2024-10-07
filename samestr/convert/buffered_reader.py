import bz2
import gzip
import tarfile


def is_gz_file(f):
    with open(f, 'rb') as _in:
        return _in.read(2) == b'\x1f\x8b'

def is_bz2_file(f):
    with open(f, 'rb') as _in:
        return _in.read(3) == b'BZh'
    
def determine_open_function(f):
    if is_gz_file(f):
        return gzip.open, "rt"
    if is_bz2_file(f):
        return bz2.BZ2File, "rb"
    return open, "rt"

def get_lines_from_chunks(stream, bufsize=800000000, binary_mode=True):
    tail = ""
    while 1:
        chunk = stream.read(bufsize)
        if binary_mode:
            chunk = chunk.decode()
        chunk = "".join((tail, chunk))            
        if not chunk:
            break
        chunk = chunk.split("\n")
        *chunk, tail = chunk
        for line in chunk:
            yield line
    if tail:
        yield tail
            
def stream_file(f):
    open_f, mode = determine_open_function(f)
    binary_mode = "t" not in mode
    with open_f(f, mode) as _in:
        yield from get_lines_from_chunks(_in, binary_mode=binary_mode)

