nextflow.enable.dsl = 2

process resource_probe {
    cpus 1
    memory '256 MB'
    time '1m'

    output:
    stdout

    script:
    """
    python3 -c 'import json; from pathlib import Path; p=Path("/sys/fs/cgroup"); print(json.dumps({"cpu_max":(p/"cpu.max").read_text().strip(),"memory_max":(p/"memory.max").read_text().strip()}))'
    """
}

workflow {
    resource_probe().view()
}
