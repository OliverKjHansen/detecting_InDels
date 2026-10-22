"""
Text files the pipeline may write gzipped: the big tables and the FASTQs, whose
names in config.yaml's `layout` end in .gz.
"""
import gzip


def open_text(path, mode="r", gzipped=None):
    """
    A text file, opened through gzip if its name ends in .gz, or if `gzipped`
    says so: a file written under a temporary name, such as x.tsv.gz.tmp,
    takes its compression from the name it will have.
    """
    path = str(path)
    if gzipped is None:
        gzipped = path.endswith(".gz")
    if not gzipped:
        return open(path, mode)
    if "r" in mode:
        return gzip.open(path, "rt")
    return gzip.open(path, mode.replace("t", "") + "t", compresslevel=6)
