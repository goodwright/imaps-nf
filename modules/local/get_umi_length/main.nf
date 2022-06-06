process GET_UMI_LENGTH {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/pysam:0.19.0--py39h5030a8b_0"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), stdout, emit: length

    script:
    def umi_separator = task.ext.args ?: ":"
    """
    #!/usr/bin/env python3

    import pysam
    import sys

    file_path = "$bam"
    umi_separator = "$umi_separator"

    bam_file = pysam.AlignmentFile(file_path, "rb")

    max_umi_len = 0
    for read in bam_file.fetch():
        if umi_separator in read.query_name:
            max_umi_len = max(max_umi_len, len(read.query_name.split(umi_separator)[-1]))

    sys.stdout.write(str(max_umi_len))
    """
}
